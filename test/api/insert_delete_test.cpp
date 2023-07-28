#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;

void initialization_args(estimate_ibf_arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
}

TEST(delete, no_given_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs_delete{0,0};
    estimate_ibf_arguments ibf_args_delete{};
    minimiser_arguments minimiser_args_delete{};
    initialization_args(ibf_args_delete);
    ibf_args_delete.compressed = false;
    ibf_args_delete.number_expression_thresholds = 2;
    minimiser_args_delete.experiment_names = false;
    ibf_args_delete.path_out = tmp_dir/"IBF_delete_Exp_";
    std::vector<std::filesystem::path> sequence_files_delete = {std::string(DATA_INPUT_DIR) + "mini_example.fasta",std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf(sequence_files_delete, ibf_args_delete, minimiser_args_delete, fpr, cutoffs_delete);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{};
    load_ibf(ibf, tmp_dir/"IBF_delete_Exp_IBF_Level_0");
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_0{seqan3::bin_count{2u},
                                         seqan3::bin_size{ibf.bin_size()},
                                         seqan3::hash_function_count{1u}};
    load_ibf(ibf, tmp_dir/"IBF_delete_Exp_IBF_Level_1");
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_1{seqan3::bin_count{2u},
                                          seqan3::bin_size{ibf.bin_size()},
                                          seqan3::hash_function_count{1u}};

    delete_bin({0,1}, ibf_args_delete,  tmp_dir/"IBF_delete_Exp_", true);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_delete{};

    load_ibf(ibf_delete, tmp_dir/"IBF_delete_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf_0 == ibf_delete));

    load_ibf(ibf_delete, tmp_dir/"IBF_delete_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf_1 == ibf_delete));

    std::filesystem::remove(tmp_dir/"IBF_delete_Exp_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_delete_Exp_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBF_delete_Exp_IBF_Deleted");
}

// Reads the level file ibf creates
template<typename float_or_int>
void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path filename)
{
    std::ifstream fin;
    fin.open(filename);
    auto stream_view = seqan3::detail::istreambuf(fin);
    auto stream_it = std::ranges::begin(stream_view);
    int j{0};
    std::vector<float_or_int> empty_vector{};

    std::string buffer{};

    // Read line = expression levels
    do
    {
        if (j == expressions.size())
            expressions.push_back(empty_vector);
        std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_char<' '>),
                                        std::back_inserter(buffer));
        if constexpr(std::same_as<uint16_t, float_or_int>)
            expressions[j].push_back((uint16_t)  std::stoi(buffer));
        else
            expressions[j].push_back((double)  std::stod(buffer));
        buffer.clear();
        if(*stream_it != '/')
            ++stream_it;

        if (*stream_it == '\n')
        {
            ++stream_it;
            j++;
        }
    } while (*stream_it != '/');
    ++stream_it;

    fin.close();
}

TEST(insert, ibf)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.compressed = false;
    ibf_args.path_out = tmp_dir/"IBF_True_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{0,0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.compressed = false;
    ibf_args_insert.path_out = tmp_dir/"IBF_True_Exp_";
    ibf_args_insert.expression_thresholds = {1, 2};
    minimiser_args_insert.experiment_names = false;
    ibf_args_insert.path_out = tmp_dir/"IBF_Insert_Exp_";
    std::vector<std::filesystem::path> sequence_files_insert = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    ibf(sequence_files_insert, ibf_args_insert, minimiser_args_insert, fpr, cutoffs_insert);

    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", tmp_dir/"IBF_Insert_Exp_", false);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
}

TEST(insert, ibf_no_given_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.compressed = false;
    ibf_args.path_out = tmp_dir/"IBF_True_Exp_";
    ibf_args.number_expression_thresholds = 2;
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05};

    std::vector<uint16_t> expected{};
    std::vector<uint8_t> cutoffs{0,0};

    std::vector<uint16_t> medians = ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.compressed = false;
    ibf_args_insert.number_expression_thresholds = 2;
    minimiser_args_insert.experiment_names = false;
    ibf_args_insert.path_out = tmp_dir/"IBF_Insert_Exp_";
    std::vector<std::filesystem::path> sequence_files_insert = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<uint16_t> medians_insert = ibf(sequence_files_insert, ibf_args_insert, minimiser_args_insert, fpr, cutoffs_insert);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", tmp_dir/"IBF_Insert_Exp_", true);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_Level_0");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, tmp_dir/"IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, tmp_dir/"IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf,expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Level_1");
}

TEST(insert, ibf_delete)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.compressed = false;
    ibf_args.path_out = tmp_dir/"IBF_True_Exp_";
    ibf_args.expression_thresholds = {1, 2};
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example2.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05,0.05};
    std::vector<uint8_t> cutoffs{0,0,0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.compressed = false;
    ibf_args_insert.path_out = tmp_dir/"IBF_Insert_Exp_";
    ibf_args_insert.expression_thresholds = {1, 2};
    minimiser_args_insert.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files_insert = {std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    std::vector<std::filesystem::path> sequence_files_test = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example2.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(sequence_files_test, ibf_args_insert, minimiser_args, fpr, cutoffs);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, false);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", tmp_dir/"IBF_Insert_Exp_", false);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_1");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_2");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Deleted");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
}


TEST(insert, ibf_delete_no_given_threshold)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    initialization_args(ibf_args);
    ibf_args.compressed = false;
    ibf_args.path_out = tmp_dir/"IBF_True_Exp_";
    ibf_args.number_expression_thresholds = 2;
    minimiser_args.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example2.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<double> fpr = {0.05,0.05};
    std::vector<uint8_t> cutoffs{0,0,0};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    std::vector<uint8_t> cutoffs_insert{0};
    estimate_ibf_arguments ibf_args_insert{};
    minimiser_arguments minimiser_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.compressed = false;
    ibf_args_insert.path_out = tmp_dir/"IBF_Insert_Exp_";
    ibf_args_insert.number_expression_thresholds = 2;
    minimiser_args_insert.experiment_names = false;
    std::vector<std::filesystem::path> sequence_files_insert = {std::string(DATA_INPUT_DIR) + "mini_example2.fasta"};
    std::vector<std::filesystem::path> sequence_files_test = {std::string(DATA_INPUT_DIR) + "mini_example.fasta", std::string(DATA_INPUT_DIR) + "mini_example2.fasta", std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    ibf(sequence_files_test, ibf_args_insert, minimiser_args, fpr, cutoffs);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, true);
    insert(sequence_files_insert, ibf_args_insert, minimiser_args_insert, cutoffs_insert, "", tmp_dir/"IBF_Insert_Exp_", true);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_Level_0");

    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, tmp_dir/"IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, tmp_dir/"IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, tmp_dir/"IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf,expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBF_True_Exp_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Deleted");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBF_Insert_Exp_IBF_FPRs.fprs");
}

TEST(insert, ibfmin)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Given_";
    ibf_args.compressed = false;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.expression_thresholds = {1, 2};
    ibf_args_insert.path_out = tmp_dir/"IBFMIN_Insert_Given_";
    ibf_args_insert.compressed = false;
    std::vector<std::filesystem::path> minimiser_file_insert = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file_insert, ibf_args_insert, fpr);
    insert(minimiser_file_insert, ibf_args_insert, "",  tmp_dir/"IBFMIN_Insert_Given_", false);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_2");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Test_Given_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_2");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
}

TEST(insert, ibfmin_delete)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.expression_thresholds = {1, 2};
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Given_";
    ibf_args.compressed = false;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.expression_thresholds = {1, 2};
    ibf_args_insert.path_out = tmp_dir/"IBFMIN_Insert_Given_";
    ibf_args_insert.compressed = false;
    std::vector<std::filesystem::path> minimiser_file_insert = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file, ibf_args_insert, fpr);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, false);
    insert(minimiser_file_insert, ibf_args_insert, "",  tmp_dir/"IBFMIN_Insert_Given_", false);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_2");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Test_Given_IBF_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_2");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Deleted");
}

TEST(insert, ibfmin_no_given_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Given_";
    ibf_args.compressed = false;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.number_expression_thresholds = 2;
    ibf_args_insert.path_out = tmp_dir/"IBFMIN_Insert_Given_";
    ibf_args_insert.compressed = false;
    fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file_insert = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file_insert, ibf_args_insert, fpr);
    insert(minimiser_file_insert, ibf_args_insert, "",  tmp_dir/"IBFMIN_Insert_Given_", true);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_Level_0");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_IBF_Level_1");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, tmp_dir/"IBFMIN_Test_Given_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf,expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_IBF_Levels.levels");
}

TEST(insert, delete_ibfmin_no_given_thresholds)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    initialization_args(ibf_args);
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05, 0.05};
    ibf_args.path_out = tmp_dir/"IBFMIN_Test_Given_Del_";
    ibf_args.compressed = false;
    std::vector<std::filesystem::path> minimiser_file = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser", std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};

    ibf(minimiser_file, ibf_args, fpr);

    estimate_ibf_arguments ibf_args_insert{};
    initialization_args(ibf_args_insert);
    ibf_args_insert.number_expression_thresholds = 2;
    ibf_args_insert.path_out = tmp_dir/"IBFMIN_Insert_Given_Del_";
    ibf_args_insert.compressed = false;
    fpr = {0.05};
    std::vector<std::filesystem::path> minimiser_file_insert = {std::string(DATA_INPUT_DIR) + "mini_example.minimiser"};
    ibf(minimiser_file, ibf_args_insert, fpr);
    delete_bin({1}, ibf_args_insert, ibf_args_insert.path_out, true);
    insert(minimiser_file_insert, ibf_args_insert, "",  tmp_dir/"IBFMIN_Insert_Given_Del_", true);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_Del_IBF_Level_0");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, tmp_dir/"IBFMIN_Test_Given_Del_IBF_Level_1");
    load_ibf(ibf_insert, tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, tmp_dir/"IBFMIN_Test_Given_Del_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf,expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, tmp_dir/"IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, tmp_dir/"IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf,fpr_insert);

    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_Del_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_Del_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_Del_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_Del_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBFMIN_Test_Given_Del_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Data");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"IBFMIN_Insert_Given_Del_IBF_Deleted");
}
