/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/rare-b-decays/b-to-k-ll-base.hh>
#include <eos/rare-b-decays/b-to-k-ll-bfs2004.hh>
#include <eos/rare-b-decays/b-to-k-ll-gp2004.hh>
#include <eos/rare-b-decays/b-to-k-ll-gvdv2020.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>


namespace eos
{
    using std::abs;
    using std::norm;
    using std::sqrt;

    struct BToKDilepton::AngularCoefficients
    {
        double a_l, b_l, c_l;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 3> & a) :
            a_l(a[0]),
            b_l(a[1]),
            c_l(a[2])
        {
        }
    };

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K} \ell^+ \ell^-@f$.
     */
    template <>
    struct Implementation<BToKDilepton>
    {
        std::shared_ptr<BToKDilepton::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        SwitchOption opt_l;

        UsedParameter hbar;
        UsedParameter m_B;
        UsedParameter m_K;
        UsedParameter m_l;
        UsedParameter tau;
        UsedParameter mu;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WilsonScan"), p, o)),
            opt_l(o, "l", { "e", "mu", "tau" }, "mu"),
            hbar(p["QM::hbar"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_K(p["mass::K_" + o.get("q", "d")], u),
            m_l(p["mass::" + opt_l.value()], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            mu(p["sb" + opt_l.value() + opt_l.value() + "::mu"], u)
        {
            std::string tag = o.get("tag", "");

            if ("BFS2004" == tag)
            {
                amplitude_generator.reset(new BToKDileptonAmplitudes<tag::BFS2004>(p, o));
            }
            else if ("GP2004" == tag)
            {
                amplitude_generator.reset(new BToKDileptonAmplitudes<tag::GP2004>(p, o));
            }
            else if ("GvDV2020" == tag)
            {
                amplitude_generator.reset(new BToKDileptonAmplitudes<tag::GvDV2020>(p, o));
            }
            else
            {
                throw InternalError("BToKDilepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*amplitude_generator);
        }

        ~Implementation()
        {
        }

        inline std::array<double, 3> angular_coefficients_array(const BToKDilepton::Amplitudes & A, const double & s) const
        {
            // cf. [BHP2007], Eq. (4.2) - (4.4)
            std::array<double, 3> result;

            // a_l
            result[0] = amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * norm(A.F_S) + norm(A.F_P))
                + 0.25 * amplitude_generator->lambda(s) * (norm(A.F_A) + norm(A.F_V))
                + 2.0 * m_l * (m_B() * m_B() - m_K() * m_K() + s) * std::real(A.F_P * std::conj(A.F_A))
                + 4.0 * m_l * m_l * m_B() * m_B() * norm(A.F_A)
                );

            // b_l
            result[1] = 2.0 * amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * std::real(A.F_S * std::conj(A.F_T))
                + std::real(A.F_P * std::conj(A.F_T5)))
                + m_l * (sqrt(amplitude_generator->lambda(s)) * beta_l(s) * std::real(A.F_S * std::conj(A.F_V))
                + (m_B() * m_B() - m_K() * m_K() + s) * std::real(A.F_T5 * std::conj(A.F_A)))
                );

            // c_l
            result[2] = amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * norm(A.F_T) + norm(A.F_T5))
                - 0.25 * amplitude_generator->lambda(s) * power_of<2>(beta_l(s)) * (norm(A.F_A) + norm(A.F_V))
                + 2.0 * m_l * sqrt(amplitude_generator->lambda(s)) * beta_l(s) * std::real(A.F_T * std::conj(A.F_V))
                );

            return result;
        }

        inline std::array<double, 3> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitude_generator->amplitudes(s), s);
        }

        inline BToKDilepton::AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return BToKDilepton::AngularCoefficients(differential_angular_coefficients_array(s));
        }

        // cf. [BHP2007], Eq. (4.8)
        inline double unnormalized_decay_width(const BToKDilepton::AngularCoefficients & a) const
        {
            return 2.0 * (a.a_l + a.c_l / 3.0);
        }

        inline double differential_branching_ratio(const BToKDilepton::AngularCoefficients & a) const
        {
            return unnormalized_decay_width(a) * tau() / hbar();
        }

        // cf. [BHP2007], Eq. (4.9)
        inline double differential_flat_term_numerator(const BToKDilepton::AngularCoefficients & a) const
        {
            return 2.0 * (a.a_l + a.c_l);
        }

        inline double differential_forward_backward_asymmetry_numerator(const BToKDilepton::AngularCoefficients & a) const
        {
            return a.b_l;
        }

        BToKDilepton::AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 3> (const double &)> integrand =
                    std::bind(&Implementation<BToKDilepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 3> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return BToKDilepton::AngularCoefficients(integrated_angular_coefficients_array);
        }

        inline double beta_l(const double & s) const
        {
            return sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

    };

    BToKDilepton::BToKDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKDilepton>(new Implementation<BToKDilepton>(parameters, options, *this))
    {
    }

    BToKDilepton::~BToKDilepton()
    {
    }

    double
    BToKDilepton::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(_imp->differential_angular_coefficients(s));
    }

    double
    BToKDilepton::differential_flat_term(const double & s) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        return _imp->differential_flat_term_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    BToKDilepton::differential_forward_backward_asymmetry(const double & s) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        return _imp->differential_forward_backward_asymmetry_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    // double
    // BToKDilepton::differential_ratio_muons_electrons(const double & s) const
    // {
    //     double br_electrons;
    //     {
    //         Save<Parameter, double> save_m_l(_imp->m_l, _imp->p["mass::e"]());
    //         Save<std::string> save_lepton_flavour(_imp->amplitude_generator->lepton_flavour, "e");
    //         br_electrons = BToKDilepton::differential_branching_ratio(s);
    //     }

    //     double br_muons;
    //     {
    //         Save<Parameter, double> save_m_l(_imp->m_l, _imp->p["mass::mu"]());
    //         Save<std::string> save_lepton_flavour(_imp->amplitude_generator->lepton_flavour, "mu");
    //         br_muons = BToKDilepton::differential_branching_ratio(s);
    //     }

    //     return br_muons / br_electrons;
    // }

    double
    BToKDilepton::two_differential_decay_width(const double & s, const double & c_theta_l_LHCb) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        // using the angular convention of the LHCb experiment
        const double c_theta_l = -c_theta_l_LHCb;

        // cf. [BHP2007], Eq. (4.1)
        return a.a_l + a.b_l * c_theta_l + a.c_l * c_theta_l * c_theta_l;
    }

    // Integrated Observables
    double
    BToKDilepton::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->unnormalized_decay_width(a);
    }

    double
    BToKDilepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_branching_ratio(a);
    }

    double
    BToKDilepton::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        double br = _imp->differential_branching_ratio(_imp->integrated_angular_coefficients(s_min, s_max));
        _imp->amplitude_generator->cp_conjugate = true;
        double br_bar = _imp->differential_branching_ratio(_imp->integrated_angular_coefficients(s_min, s_max));

        return (br + br_bar) / 2.0;
    }

    double
    BToKDilepton::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        auto gamma     = _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));
        _imp->amplitude_generator->cp_conjugate = true;
        auto gamma_bar = _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }

    double
    BToKDilepton::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_flat_term_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    BToKDilepton::integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        double num_integrated = _imp->differential_flat_term_numerator(_imp->integrated_angular_coefficients(s_min, s_max));
        double denom_integrated = _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        _imp->amplitude_generator->cp_conjugate = true;

        num_integrated += _imp->differential_flat_term_numerator(_imp->integrated_angular_coefficients(s_min, s_max));
        denom_integrated += _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        return num_integrated / denom_integrated;
    }

    double
    BToKDilepton::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_forward_backward_asymmetry_numerator(a) / _imp->unnormalized_decay_width(a);

    }

    double
    BToKDilepton::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        double num_integrated = _imp->differential_forward_backward_asymmetry_numerator(_imp->integrated_angular_coefficients(s_min, s_max));
        double denom_integrated = _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        _imp->amplitude_generator->cp_conjugate = true;

        num_integrated += _imp->differential_forward_backward_asymmetry_numerator(_imp->integrated_angular_coefficients(s_min, s_max));
        denom_integrated += _imp->unnormalized_decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        return num_integrated / denom_integrated;
    }

    // double
    // BToKDilepton::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    // {
    //     std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToKDilepton::differential_branching_ratio),
    //             this, std::placeholders::_1);

    //     double br_electrons;
    //     {
    //         Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
    //         Save<std::string> save_lepton_flavour(_imp->amplitude_generator->lepton_flavour, "e");
    //         // br_electrons = integrate<GSL::QNG>(integrand, s_min, s_max);
    //         br_electrons = integrate<GSL::QNG>(integrand, s_min, s_max);
    //     }

    //     double br_muons;
    //     {
    //         Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
    //         Save<std::string> save_lepton_flavour(_imp->amplitude_generator->lepton_flavour, "mu");
    //         br_muons = integrate<GSL::QNG>(integrand, s_min, s_max);
    //     }

    //     // cf. [BHP2007], Eq. (4.10), p. 6
    //     return br_muons / br_electrons;
    // }


    const std::string
    BToKDilepton::description = "\
The decay B->K l^+ l^-, with l=e,mu,tau a charged lepton.";

    const std::string
    BToKDilepton::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    BToKDilepton::kinematics_description_c_theta_l = "\
The cosine of the lepton's helicity angle theta_l in the l^+l^- rest frame using the LHCb convention.";

    /*
     * For diagnostic purposes only!
     */
    BToKDilepton::Amplitudes
    BToKDilepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }

    std::array<double, 3>
    BToKDilepton::angular_coefficients(const double & q2) const
    {
        return _imp->angular_coefficients_array(_imp->amplitude_generator->amplitudes(q2), q2);
    }
}
