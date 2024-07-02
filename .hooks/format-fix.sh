#!/bin/bash

cmake --build ${BUILD_DIR} --parallel 12 --config Debug --target format-fix
exit $?
