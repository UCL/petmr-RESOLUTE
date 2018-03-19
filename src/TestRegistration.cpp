/*
   Resolute.cpp

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

   This program tests RESOLUTE mask extraction.
 */

#include <iostream>
#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>
#include <nlohmann/json.hpp>

#include "EnvironmentInfo.h"

#include "antsRegistrationTemplateHeader.h"
#include <include/ants.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

namespace fs = boost::filesystem;
using namespace ants;

int main(int argc, char **argv)
{

std::string argList = "3 -m CC[<%%REF%%>,<%%FLOAT%%>,1,4] -i 10x5x2 -o <%%PREFIX%%> -t SyN[0.5] -r Gauss[3,0] -G";

std::vector<std::string> args;

//const boost::regex ex("--");
//boost::split(args, argList, ex, boost::token_compress_on);

boost::split_regex( args, argList, boost::regex( " " ) ) ;

int x=1;
for (auto a : args){
  LOG(INFO) << x << "\t\t" << a;
  x++;
}

ants::ANTS( args, &std::cout);

  return EXIT_SUCCESS;
}