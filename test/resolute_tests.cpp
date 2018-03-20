/*
   resolute_tests.cpp

   Author:      Benjamin A. Thomas

   Copyright 2018 Institute of Nuclear Medicine, University College London.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
 */

#include "Resolute.hpp"
#include <gtest/gtest.h>

namespace {

TEST(Resolute, R2StoHU){
   EXPECT_NEAR(330, ns::GetHUfromR2s(100.0), 10);
   EXPECT_NEAR(615, ns::GetHUfromR2s(200.0), 10);
   EXPECT_NEAR(844, ns::GetHUfromR2s(300.0), 10);
   EXPECT_NEAR(1025, ns::GetHUfromR2s(400.0), 10);
   EXPECT_NEAR(1166, ns::GetHUfromR2s(500.0), 10);
   EXPECT_NEAR(1275, ns::GetHUfromR2s(600.0), 10);
   EXPECT_NEAR(1360, ns::GetHUfromR2s(700.0), 10);
   EXPECT_NEAR(1430, ns::GetHUfromR2s(800.0), 10);
   EXPECT_NEAR(1492, ns::GetHUfromR2s(900.0), 10);
   EXPECT_NEAR(1556, ns::GetHUfromR2s(1000.0), 10);
}

TEST(Resolute, R2StoLAC){
   EXPECT_NEAR(0.115, ns::GetMU(100.0), 0.001);
   EXPECT_NEAR(0.129, ns::GetMU(200.0), 0.001);
   EXPECT_NEAR(0.141, ns::GetMU(300.0), 0.001);
   EXPECT_NEAR(0.150, ns::GetMU(400.0), 0.001);
   EXPECT_NEAR(0.158, ns::GetMU(500.0), 0.001);
   EXPECT_NEAR(0.163, ns::GetMU(600.0), 0.001);
   EXPECT_NEAR(0.167, ns::GetMU(700.0), 0.001);
   EXPECT_NEAR(0.171, ns::GetMU(800.0), 0.001);
   EXPECT_NEAR(0.174, ns::GetMU(900.0), 0.001);
   EXPECT_NEAR(0.177, ns::GetMU(1000.0), 0.001);
}

}