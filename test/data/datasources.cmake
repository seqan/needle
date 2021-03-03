cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/<file>
declare_datasource (FILE exp_01.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/exp_01.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
declare_datasource (FILE IBF_1
                    URL ${CMAKE_SOURCE_DIR}/test/data/IBF_1
                    URL_HASH SHA256=e752c7d281ffa8a5104c6c6bdfd0c4868d39cf98bb26e95c4adbfe52c9eb6784)
declare_datasource (FILE mini_example.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.fasta
                    URL_HASH SHA256=f872221632e7423b071d1aedaf4c4e9da2a659e843fcafac12cd65632b904b93)
declare_datasource (FILE mini_example.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.header
                    URL_HASH SHA256=54644375bf8e1d3d9d8ba556762afd3f88f5eba5ca4324629b79517b7aea36a1)
declare_datasource (FILE mini_example.minimiser
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.minimiser
                    URL_HASH SHA256=98927a56465368db8c8a557f0ad5b83c25f57fb1820ba997be23b69fb4fe9244)
declare_datasource (FILE mini_example2.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example2.header
                    URL_HASH SHA256=095762f665803679c6149e0dbad2462526f68c52f477a4aaf21386ad82aaea07)
declare_datasource (FILE mini_gen.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen.fasta
                    URL_HASH SHA256=6e9da2f6693938586c902f5e4445f2df1f1ac94cff8c23dea9e02b58759a8998)
declare_datasource (FILE mini_gen2.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen2.fasta
                    URL_HASH SHA256=7e21b3eff20f950e6a0fa2369a0bbbed4cc967a24083876350f839f4a232771b)
declare_datasource (FILE mini_genom.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_genom.fasta
                    URL_HASH SHA256=b53541bb438be5241880828048fd37544c6aca30d363db46d0d918a2bc531a0e)
