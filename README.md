# SeqAn3 App Template

[![Build Status](https://travis-ci.com/joergi-w/app-template.svg?branch=master)](https://travis-ci.com/joergi-w/app-template)

This is a template for app developers with SeqAn3. You can easily clone this repository and modify the existing code to your needs. It provides the elementary set-up for all SeqAn3 applications.

Instructions:
1. clone this repository: `git clone --recurse-submodules https://github.com/joergi-w/app-template.git`
2. create a build directory and visit it: `mkdir build && cd build`
3. run cmake: `cmake ../app-template`
4. run `make`
5. run the tests: `make test`
6. execute the app: `./my_app`

Planned features are:
- build the app with seqan3
- cli tests (black box)
- api tests (white box)
- tutorials, usage instructions
- api documentation
- benchmarks for time, memory and disk space consumption
- coverage tests: api
- automated CD
- automated CI
