/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Danny van Dyk
 * Copyright (c) 2021 Stefan Meiser
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

#include <eos/b-decays/b-to-d-psd.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/options-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: B_q -> D_q P, cf. [BBNS:2000A]
     */
    template <>
    struct Implementation<BToDPseudoscalar>
    {
        std::shared_ptr<Model> model;

        SwitchOption opt_q;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_B;

        UsedParameter tau_B;

        bool cp_conjugate;

        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, "q", { "d", "s" }),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_" + opt_q.value()], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            mu(p["cbsu::mu"], u)
        {
        }

        double decay_width() const
        {
            //const auto wc = model->wet_cbsu(cp_conjugate);

            // masses
            const double m_B = this-> m_B(), m_B2 = m_B * m_B;

            return power_of<2>(g_fermi * std::abs(model->ckm_cb()))
                * 1.0;
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToDPseudoscalar>::options
    {
	    { "q", { "d", "s" }, ""}
    };

    BToDPseudoscalar::BToDPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDPseudoscalar>(new Implementation<BToDPseudoscalar>(parameters, options, *this))
    {
    }

    BToDPseudoscalar::~BToDPseudoscalar()
    {
    }

    double
    BToDPseudoscalar::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    const std::set<ReferenceName>
    BToDPseudoscalar::references
    {
        "BBNS:2000A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::begin_options()
    {
        return Implementation<BToDPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToDPseudoscalar::end_options()
    {
        return Implementation<BToDPseudoscalar>::options.cend();
    }
}
