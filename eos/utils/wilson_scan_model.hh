/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013, 2015 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014, 2018 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_SRC_UTILS_WILSON_SCAN_MODEL_HH
#define EOS_GUARD_SRC_UTILS_WILSON_SCAN_MODEL_HH 1

#include <eos/utils/ckm_scan_model.hh>
#include <eos/utils/model.hh>
#include <eos/utils/standard-model.hh>

namespace eos
{
    template <typename Tag_> class WilsonScanComponent;

    template <>
    class WilsonScanComponent<components::DeltaBS1> :
        public virtual ModelComponent<components::DeltaBS1>
    {
        protected:
            /* QCD parameters */
            UsedParameter _alpha_s_Z__deltabs1;
            UsedParameter _mu_b__deltabs1;

            /* Masses */
            UsedParameter _m_Z__deltabs1;

            /* Renormalization scale */
            UsedParameter _mu__deltabs1;

            /* b->s Wilson coefficients */
            UsedParameter _c1;
            UsedParameter _c2;
            UsedParameter _c3;
            UsedParameter _c4;
            UsedParameter _c5;
            UsedParameter _c6;
            UsedParameter _re_c7,          _im_c7;
            UsedParameter _re_c7prime,     _im_c7prime;
            UsedParameter _c8;
            UsedParameter _c8prime;
            /* b->see Wilson coefficients */
            UsedParameter _e_re_c9,        _e_im_c9;
            UsedParameter _e_re_c10,       _e_im_c10;
            UsedParameter _e_re_c9prime,   _e_im_c9prime;
            UsedParameter _e_re_c10prime,  _e_im_c10prime;
            UsedParameter _e_re_cS,        _e_im_cS;
            UsedParameter _e_re_cSprime,   _e_im_cSprime;
            UsedParameter _e_re_cP,        _e_im_cP;
            UsedParameter _e_re_cPprime,   _e_im_cPprime;
            UsedParameter _e_re_cT,        _e_im_cT;
            UsedParameter _e_re_cT5,       _e_im_cT5;
            /* b->smumu Wilson coefficients */
            UsedParameter _mu_re_c9,       _mu_im_c9;
            UsedParameter _mu_re_c10,      _mu_im_c10;
            UsedParameter _mu_re_c9prime,  _mu_im_c9prime;
            UsedParameter _mu_re_c10prime, _mu_im_c10prime;
            UsedParameter _mu_re_cS,       _mu_im_cS;
            UsedParameter _mu_re_cSprime,  _mu_im_cSprime;
            UsedParameter _mu_re_cP,       _mu_im_cP;
            UsedParameter _mu_re_cPprime,  _mu_im_cPprime;
            UsedParameter _mu_re_cT,       _mu_im_cT;
            UsedParameter _mu_re_cT5,      _mu_im_cT5;

            /* b->sgamma */
            std::function<complex<double> ()> _c7;
            std::function<complex<double> ()> _c7prime;

            /* b->see */
            std::function<complex<double> ()> _e_c9;
            std::function<complex<double> ()> _e_c10;
            std::function<complex<double> ()> _e_c9prime;
            std::function<complex<double> ()> _e_c10prime;
            std::function<complex<double> ()> _e_cS;
            std::function<complex<double> ()> _e_cSprime;
            std::function<complex<double> ()> _e_cP;
            std::function<complex<double> ()> _e_cPprime;
            std::function<complex<double> ()> _e_cT;
            std::function<complex<double> ()> _e_cT5;

            /* b->smumu */
            std::function<complex<double> ()> _mu_c9;
            std::function<complex<double> ()> _mu_c10;
            std::function<complex<double> ()> _mu_c9prime;
            std::function<complex<double> ()> _mu_c10prime;
            std::function<complex<double> ()> _mu_cS;
            std::function<complex<double> ()> _mu_cSprime;
            std::function<complex<double> ()> _mu_cP;
            std::function<complex<double> ()> _mu_cPprime;
            std::function<complex<double> ()> _mu_cT;
            std::function<complex<double> ()> _mu_cT5;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /*! b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s(const double & mu, const std::string & lepton_flavor, const bool & cp_conjugate) const;
    };

    template <>
    class WilsonScanComponent<components::WET::SBSB> :
        public virtual ModelComponent<components::WET::SBSB>
    {
        protected:
            /* b->s Wilson coefficients */
            UsedParameter _re_sbsb_c1__deltab2;
            UsedParameter _im_sbsb_c1__deltab2;
            UsedParameter _re_sbsb_c2__deltab2;
            UsedParameter _im_sbsb_c2__deltab2;
            UsedParameter _re_sbsb_c3__deltab2;
            UsedParameter _im_sbsb_c3__deltab2;
            UsedParameter _re_sbsb_c4__deltab2;
            UsedParameter _im_sbsb_c4__deltab2;
            UsedParameter _re_sbsb_c5__deltab2;
            UsedParameter _im_sbsb_c5__deltab2;
            UsedParameter _re_sbsb_c1p__deltab2;
            UsedParameter _im_sbsb_c1p__deltab2;
            UsedParameter _re_sbsb_c2p__deltab2;
            UsedParameter _im_sbsb_c2p__deltab2;
            UsedParameter _re_sbsb_c3p__deltab2;
            UsedParameter _im_sbsb_c3p__deltab2;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /*! sbar b sbar b Wilson coefficients */
            virtual WilsonCoefficients<wc::SBSB> wet_sbsb() const;
    };

    template <>
    class WilsonScanComponent<components::WET::UBLNu> :
        public virtual ModelComponent<components::WET::UBLNu>
    {
        private:
            /* b->u Wilson coefficients */
            /* b->u e nu_e */
            UsedParameter _e_re_csl, _e_im_csl;
            UsedParameter _e_re_csr, _e_im_csr;
            UsedParameter _e_re_cvl, _e_im_cvl;
            UsedParameter _e_re_cvr, _e_im_cvr;
            UsedParameter _e_re_ct,  _e_im_ct;

            /* b->u mu nu_mu */
            UsedParameter _mu_re_csl, _mu_im_csl;
            UsedParameter _mu_re_csr, _mu_im_csr;
            UsedParameter _mu_re_cvl, _mu_im_cvl;
            UsedParameter _mu_re_cvr, _mu_im_cvr;
            UsedParameter _mu_re_ct,  _mu_im_ct;

            /* b->u tau nu_tau */
            UsedParameter _tau_re_csl, _tau_im_csl;
            UsedParameter _tau_re_csr, _tau_im_csr;
            UsedParameter _tau_re_cvl, _tau_im_cvl;
            UsedParameter _tau_re_cvr, _tau_im_cvr;
            UsedParameter _tau_re_ct,  _tau_im_ct;

            std::function<complex<double> ()> _e_csl;
            std::function<complex<double> ()> _e_csr;
            std::function<complex<double> ()> _e_cvl;
            std::function<complex<double> ()> _e_cvr;
            std::function<complex<double> ()> _e_ct;

            std::function<complex<double> ()> _mu_csl;
            std::function<complex<double> ()> _mu_csr;
            std::function<complex<double> ()> _mu_cvl;
            std::function<complex<double> ()> _mu_cvr;
            std::function<complex<double> ()> _mu_ct;

            std::function<complex<double> ()> _tau_csl;
            std::function<complex<double> ()> _tau_csr;
            std::function<complex<double> ()> _tau_cvl;
            std::function<complex<double> ()> _tau_cvr;
            std::function<complex<double> ()> _tau_ct;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /* b->u Wilson coefficients */
            virtual WilsonCoefficients<ChargedCurrent> wet_ublnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const;
    };

    template <>
    class WilsonScanComponent<components::WET::CBLNu> :
    public virtual ModelComponent<components::WET::CBLNu>
    {
        private:
            /* b->c Wilson coefficients */
            /* b->c e nu_e */
            UsedParameter _e_re_csl, _e_im_csl;
            UsedParameter _e_re_csr, _e_im_csr;
            UsedParameter _e_re_cvl, _e_im_cvl;
            UsedParameter _e_re_cvr, _e_im_cvr;
            UsedParameter _e_re_ct,  _e_im_ct;

            /* b->c mu nu_mu */
            UsedParameter _mu_re_csl, _mu_im_csl;
            UsedParameter _mu_re_csr, _mu_im_csr;
            UsedParameter _mu_re_cvl, _mu_im_cvl;
            UsedParameter _mu_re_cvr, _mu_im_cvr;
            UsedParameter _mu_re_ct,  _mu_im_ct;

            /* b->c tau nu_tau */
            UsedParameter _tau_re_csl, _tau_im_csl;
            UsedParameter _tau_re_csr, _tau_im_csr;
            UsedParameter _tau_re_cvl, _tau_im_cvl;
            UsedParameter _tau_re_cvr, _tau_im_cvr;
            UsedParameter _tau_re_ct,  _tau_im_ct;

            std::function<complex<double> ()> _e_csl;
            std::function<complex<double> ()> _e_csr;
            std::function<complex<double> ()> _e_cvl;
            std::function<complex<double> ()> _e_cvr;
            std::function<complex<double> ()> _e_ct;

            std::function<complex<double> ()> _mu_csl;
            std::function<complex<double> ()> _mu_csr;
            std::function<complex<double> ()> _mu_cvl;
            std::function<complex<double> ()> _mu_cvr;
            std::function<complex<double> ()> _mu_ct;

            std::function<complex<double> ()> _tau_csl;
            std::function<complex<double> ()> _tau_csr;
            std::function<complex<double> ()> _tau_cvl;
            std::function<complex<double> ()> _tau_cvr;
            std::function<complex<double> ()> _tau_ct;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /* b->c Wilson coefficients */
            virtual WilsonCoefficients<ChargedCurrent> wet_cblnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate) const;
    };

    template <>
    class WilsonScanComponent<components::WET::SBNuNu> :
    public virtual ModelComponent<components::WET::SBNuNu>
    {
        private:
            /* sbnunu Wilson coefficients */
            UsedParameter _re_cl, _im_cl;
            UsedParameter _re_cr, _im_cr;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /* sbnunu Wilson coefficients */
            virtual WilsonCoefficients<wc::SBNuNu> wet_sbnunu(const bool & cp_conjugate) const;
    };

    template <>
    class WilsonScanComponent<components::CBSU> :
        public virtual ModelComponent<components::CBSU>
    {
        protected:
            /* b->s c u Wilson coefficients */
            UsedParameter _re_cbsu_c1__CBSU;
            UsedParameter _im_cbsu_c1__CBSU;
            UsedParameter _re_cbsu_c2__CBSU;
            UsedParameter _im_cbsu_c2__CBSU;
            UsedParameter _re_cbsu_c3__CBSU;
            UsedParameter _im_cbsu_c3__CBSU;
            UsedParameter _re_cbsu_c4__CBSU;
            UsedParameter _im_cbsu_c4__CBSU;
            UsedParameter _re_cbsu_c5__CBSU;
            UsedParameter _im_cbsu_c5__CBSU;
            UsedParameter _re_cbsu_c6__CBSU;
            UsedParameter _im_cbsu_c6__CBSU;
            UsedParameter _re_cbsu_c7__CBSU;
            UsedParameter _im_cbsu_c7__CBSU;
            UsedParameter _re_cbsu_c8__CBSU;
            UsedParameter _im_cbsu_c8__CBSU;
            UsedParameter _re_cbsu_c9__CBSU;
            UsedParameter _im_cbsu_c9__CBSU;
            UsedParameter _re_cbsu_c10__CBSU;
            UsedParameter _im_cbsu_c10__CBSU;
            UsedParameter _re_cbsu_c1p__CBSU;
            UsedParameter _im_cbsu_c1p__CBSU;
            UsedParameter _re_cbsu_c2p__CBSU;
            UsedParameter _im_cbsu_c2p__CBSU;
            UsedParameter _re_cbsu_c3p__CBSU;
            UsedParameter _im_cbsu_c3p__CBSU;
            UsedParameter _re_cbsu_c4p__CBSU;
            UsedParameter _im_cbsu_c4p__CBSU;
            UsedParameter _re_cbsu_c5p__CBSU;
            UsedParameter _im_cbsu_c5p__CBSU;
            UsedParameter _re_cbsu_c6p__CBSU;
            UsedParameter _im_cbsu_c6p__CBSU;
            UsedParameter _re_cbsu_c7p__CBSU;
            UsedParameter _im_cbsu_c7p__CBSU;
            UsedParameter _re_cbsu_c8p__CBSU;
            UsedParameter _im_cbsu_c8p__CBSU;
            UsedParameter _re_cbsu_c9p__CBSU;
            UsedParameter _im_cbsu_c9p__CBSU;
            UsedParameter _re_cbsu_c10p__CBSU;
            UsedParameter _im_cbsu_c10p__CBSU;

        public:
            WilsonScanComponent(const Parameters &, const Options &, ParameterUser &);

            /*! cbar b sbar u Wilson coefficients */
            virtual WilsonCoefficients<wc::CBSU> wilson_coefficients_cbsu(const double & mu) const;
    };

    /*!
     * A model with all possible operators; their Wilson coefficients
     * are allowed to have arbitrary values.
     */
    class WilsonScanModel :
        public Model,
        public CKMScanComponent,
        public SMComponent<components::QCD>,
        public WilsonScanComponent<components::WET::SBSB>,
        public WilsonScanComponent<components::DeltaBS1>,
        public WilsonScanComponent<components::WET::UBLNu>,
        public WilsonScanComponent<components::WET::CBLNu>,
        public WilsonScanComponent<components::WET::SBNuNu>,
        public WilsonScanComponent<components::WET::CBSU>
    {
        public:
            WilsonScanModel(const Parameters &, const Options &);
            virtual ~WilsonScanModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };

    class ConstrainedWilsonScanComponent :
        public WilsonScanComponent<components::DeltaBS1>
    {
        public:
            ConstrainedWilsonScanComponent(const Parameters &, const Options &, ParameterUser &);
    };

    /*!
     * Special case of @see WilsonScanModel with C_S = - C_P, C'_S = C'_P, and C_T = C_T5 = 0.
     *
     * As shown in arXiv:1407.7044 eq. (8), the Wilson coefficients are not
     * independent if new physics is well above the electro-weak scale,
     * respects the SM gauge symmetry, and only dim. 6 operators contribute.
     */
    class ConstrainedWilsonScanModel :
        public Model,
        public CKMScanComponent,
        public SMComponent<components::QCD>,
        public WilsonScanComponent<components::WET::SBSB>,
        public ConstrainedWilsonScanComponent,
        public WilsonScanComponent<components::WET::UBLNu>,
        public WilsonScanComponent<components::WET::CBLNu>,
        public WilsonScanComponent<components::WET::SBNuNu>,
        public WilsonScanComponent<components::WET::CBSU>
    {
        public:
            ConstrainedWilsonScanModel(const Parameters &, const Options &);
            virtual ~ConstrainedWilsonScanModel();

            static std::shared_ptr<Model> make(const Parameters &, const Options &);
    };
}

#endif
