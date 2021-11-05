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

#include <test/test.hh>
#include <eos/observable.hh>
#include <eos/b-decays/b-to-d-psd.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>


using namespace test;
using namespace eos;

class BToDPseudoscalarTest :
    public TestCase
{
    public:
        BToDPseudoscalarTest() :
            TestCase("b_to_d_psd_test")
        {
        }

        virtual void run() const
        {
            {
                Parameters p = Parameters::Defaults();

                Options oo
                {
                    { "model", "CKM" },
                    { "q",     "s"   }
                };

                BToDPseudoscalar d(p, oo);

                const double eps = 1e-12;

                TEST_CHECK_NEARLY_EQUAL(
                        1.000,
                        d.branching_ratio(),
                        eps);
            }
        }
} b_to_d_psd_test;
