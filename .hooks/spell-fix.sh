#!/bin/bash

cmake --build ${BUILD_DIR} --parallel 12 --config Debug --target spell-fix
exit $?
