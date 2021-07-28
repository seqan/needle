cmake_minimum_required (VERSION 3.9)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/<file>
declare_datasource (FILE exp_01.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/exp_01.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
declare_datasource (FILE IBF_1
                    URL ${CMAKE_SOURCE_DIR}/test/data/IBF_1
                    URL_HASH SHA256=b0f23b093f774e29f75d03c7ebb94c5e809ee2f640d6b6d0e93667ad251ad480)
declare_datasource (FILE mini_example.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.fasta
                    URL_HASH SHA256=f872221632e7423b071d1aedaf4c4e9da2a659e843fcafac12cd65632b904b93)
declare_datasource (FILE mini_example.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.header
                    URL_HASH SHA256=f31019786eff444ba5fd1f8a04df070b4dda1269b19aa8ad285e1cc61e2b25ce)
declare_datasource (FILE mini_example.minimiser
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.minimiser
                    URL_HASH SHA256=7848a3c0eb438cd940ad21127e836c9b82a97ab1e92ed1e0102f392ab60d5f1b)
declare_datasource (FILE mini_example2.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example2.header
                    URL_HASH SHA256=9df57d832dea68fe2ef19d783337f9193521b3ad3664b6431d79cd807a362e9f)
declare_datasource (FILE mini_gen.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen.fasta
                    URL_HASH SHA256=6e9da2f6693938586c902f5e4445f2df1f1ac94cff8c23dea9e02b58759a8998)
declare_datasource (FILE mini_gen2.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen2.fasta
                    URL_HASH SHA256=7e21b3eff20f950e6a0fa2369a0bbbed4cc967a24083876350f839f4a232771b)
