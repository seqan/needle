metagraph="Set path to metagraph executable"
out_dir="Set out directory"

mkdir $out_dir/clean_anno
/usr/bin/time -v -o jointbuild.time  $metagraph build -p 4 -k 21 --count-kmers --count-width 32 --outfile-base $out_dir/graph single_dbgs/clean/SRR*.fasta.gz &> jointbuilt.log
/usr/bin/time -v -o annotate_joint.time  $metagraph annotate --separately -p 4 -i $out_dir/graph.dbg --anno-filename --count-kmers --count-width 32 -o $out_dir/clean_anno single_dbgs/clean/SRR*.fasta.gz &> anno.log

mkdir $out_dir/rd0
mkdir $out_dir/rd1
mkdir $out_dir/rd2
/usr/bin/time -v -o rd0.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 0 --mem-cap-gb 500 -o $out_dir/rd0/rd -i $out_dir/graph.dbg -p 4 --disk-swap "" $out_dir/clean_anno/*.column.annodbg &>rd0.log
/usr/bin/time -v -o rd1.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 1 --mem-cap-gb 500 -o $out_dir/rd1/rd -i $out_dir/graph.dbg -p 4 --disk-swap "" $out_dir/clean_anno/*.column.annodbg &> rd1.log
/usr/bin/time -v -o rd2.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 2 --mem-cap-gb 500 -o $out_dir/rd2/rd -i $out_dir/graph.dbg -p 4 --disk-swap "" $out_dir/clean_anno/*.column.annodbg &> rd2.log

/usr/bin/time -v -o anno.time $metagraph transform_anno --anno-type row_diff_int_brwt --greedy --fast --subsample 1000000 -i $out_dir/graph.dbg -o $out_dir/annotation_final $out_dir/rd2/*.column.annodbg -p 4 --parallel-nodes 10 &> last_anno.log
/usr/bin/time -v -o relax.time $metagraph relax_brwt -v --relax-arity 32 -p 4 -o $out_dir/annotation_final_relaxed $out_dir/annotation_final.row_diff_int_brwt.annodbg &> relax.log

# Note data/*.fa are the files from here: https://github.com/MitraDarja/analysis_needle/tree/main
/usr/bin/time -v -o query_1.time  $metagraph query -p 4 --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_1.fa
/usr/bin/time -v -o query_100.time  $metagraph query -p 4  --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_100.fa
/usr/bin/time -v -o query_1000.time  $metagraph query -p 4 --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_1000.fa

# Query with one thead
/usr/bin/time -v -o query_1_1.time  $metagraph query -p 1 --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_1.fa
/usr/bin/time -v -o query_100_1.time  $metagraph query -p 1  --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_100.fa
/usr/bin/time -v -o query_1000_1.time  $metagraph query -p 1 --query-counts -i $out_dir/graph.dbg -a $out_dir/annotation_final_relaxed.row_diff_int_brwt.annodbg data/query_1000.fa


# Construction for smooth option

mkdir $out_dir/clean_anno_smooth
/usr/bin/time -v -o jointbuild_smooth.time  $metagraph build -p 4 -k 21 --count-kmers --count-width 32 --outfile-base $out_dir/graph_smooth single_dbgs/clean_smooth/SRR*.fasta.gz &> jointbuilt.log
/usr/bin/time -v -o annotate_joint_smooth.time  $metagraph annotate --separately -p 4 -i $out_dir/graph.dbg --anno-filename --count-kmers --count-width 32 -o $out_dir/clean_anno_smooth single_dbgs/clean_smooth/SRR*.fasta.gz &> anno.log

mkdir $out_dir/smooth_rd0
mkdir $out_dir/smooth_rd1
mkdir $out_dir/smooth_rd2
/usr/bin/time -v -o rd0_smooth.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 0 --mem-cap-gb 500 -o $out_dir/smooth_rd0/rd -i $out_dir/graph_smooth.dbg -p 4 --disk-swap "" $out_dir/clean_anno_smooth/*.column.annodbg &>rd0.log
/usr/bin/time -v -o rd1_smooth.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 1 --mem-cap-gb 500 -o $out_dir/smooth_rd1/rd -i $out_dir/graph_smooth.dbg -p 4 --disk-swap "" $out_dir/clean_anno_smooth/*.column.annodbg &> rd1.log
/usr/bin/time -v -o rd2_smooth.time  $metagraph transform_anno -v --anno-type row_diff --count-kmers --row-diff-stage 2 --mem-cap-gb 500 -o $out_dir/smooth_rd2/rd -i $out_dir/graph_smooth.dbg -p 4 --disk-swap "" $out_dir/clean_anno_smooth/*.column.annodbg &> rd2.log

/usr/bin/time -v -o anno_smooth.time $metagraph transform_anno --anno-type row_diff_int_brwt --greedy --fast --subsample 1000000 -i $out_dir/graph_smooth.dbg -o $out_dir/annotation_final_smooth $out_dir/rd2/*.column.annodbg -p 4 --parallel-nodes 10 &> last_anno.log
/usr/bin/time -v -o relax_smooth.time $metagraph relax_brwt -v --relax-arity 32 -p 4 -o $out_dir/annotation_final_relaxed_smooth $out_dir/annotation_final_smooth.row_diff_int_brwt.annodbg &> relax.log

# Note data/*.fa are the files from here: https://github.com/MitraDarja/analysis_needle/tree/main
/usr/bin/time -v -o query_1_smooth.time  $metagraph query -p 4 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_1.fa
/usr/bin/time -v -o query_100_smooth.time  $metagraph query -p 4 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_100.fa
/usr/bin/time -v -o query_1000_smooth.time  $metagraph query -p 4 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_1000.fa

# Query with one thead
/usr/bin/time -v -o query_1_smooth_1.time  $metagraph query -p 1 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_1.fa
/usr/bin/time -v -o query_100_smooth_1.time  $metagraph query -p 1 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_100.fa
/usr/bin/time -v -o query_1000_smooth_1.time  $metagraph query -p 1 --query-counts -i $out_dir/graph_smooth.dbg -a $out_dir/annotation_final_relaxed_smooth.row_diff_int_brwt.annodbg data/query_1000.fa
