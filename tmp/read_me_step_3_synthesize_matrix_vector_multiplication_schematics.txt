Calling convention

synthesator [<system name>='matrix' [<GF2>='0' [<pipelinedelay>='0' <timescale>='1']]]
or
synthesator [<matrix file name (default matrix)> [<GF2(default 0)> [<pipelinedelay(default 0> <timescale (default 1)>]]]

here:

<system name> - char* path and name of a *.txt file with original matrix without .txt extention

<GF2> - int, 
	= 0 - arithmetic matrix (*,+)
	< 0 - geometric matrix (^(0,1),*) - needs hand corrections if kernel blocking is used
	> 0 - logical matrix (&,xor) - - needs hand corrections if kernel blocking is used

<piprlinedalay> - >=0, delays in sum or mult operators if GF2<=0

<timescale> - int,
	<= 0 - all balansing delays are 0, idempotent operation, inputs must be separated.
	> 0	- serial input, time depenency not destorted


examples

arithmetic, pipeline delay = 2, serial input:

..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\B_comp 0 2 1
..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\Bt_comp 0 2 1
..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\Pr_comp 0 2 1


geometric, pipeline delay = 3, serial input:

..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\S1_comp -1 3 1


logical with blocking = 8, pipeline delay = 0 (always), serial input

..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\GenMat_comp 8 0 1
..\step_2_compile_code_for_matrix_vector_multiplier\synthesator_suboptimal_web_c_code\x64\Release\synthesator.exe ..\step_1_generate_compact_matrices\H_comp 8 0 1 
