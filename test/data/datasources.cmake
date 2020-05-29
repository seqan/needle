cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/in.fastq
declare_datasource (FILE exp_01.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/exp_01.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
declare_datasource (FILE IBF_1.000000
                    URL ${CMAKE_SOURCE_DIR}/test/data/IBF_1.000000
                    URL_HASH SHA256=d43e898adada7a656cbfdcc8b88c024f8e1a4ad6ee31d355ddac11e04c896aac)
declare_datasource (FILE mini_example.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.fasta
                    URL_HASH SHA256=f872221632e7423b071d1aedaf4c4e9da2a659e843fcafac12cd65632b904b93)
declare_datasource (FILE mini_example.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.header
                    URL_HASH SHA256=685797f30f29436e74b4e373329dba61afcb88b6d73e490754aba8b3f3d78002)
declare_datasource (FILE mini_example.minimiser
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.minimiser
                    URL_HASH SHA256=98927a56465368db8c8a557f0ad5b83c25f57fb1820ba997be23b69fb4fe9244)
declare_datasource (FILE mini_example2.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example2.header
                    URL_HASH SHA256=9cd8e62254bfd8313930845d52ca95fbaaf2cdb70b45367b6d04ce9b25382403)
declare_datasource (FILE mini_gen.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen.fasta
                    URL_HASH SHA256=6e9da2f6693938586c902f5e4445f2df1f1ac94cff8c23dea9e02b58759a8998)
declare_datasource (FILE mini_gen2.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen2.fasta
                    URL_HASH SHA256=7e21b3eff20f950e6a0fa2369a0bbbed4cc967a24083876350f839f4a232771b)
declare_datasource (FILE mini_genom.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_genom.fasta
                    URL_HASH SHA256=b53541bb438be5241880828048fd37544c6aca30d363db46d0d918a2bc531a0e)
