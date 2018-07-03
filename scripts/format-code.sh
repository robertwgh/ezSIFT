#!/bin/bash
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
clang-format -i ${dir}/../**/*.cpp ${dir}/../**/*.h


