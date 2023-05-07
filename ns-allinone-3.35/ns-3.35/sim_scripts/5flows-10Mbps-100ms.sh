#!/bin/bash
cd ..

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=0 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=0 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=0 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=0 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=346 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=346 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=346 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=346 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=519 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=519 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=519 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=519 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=692 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=692 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=692 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=692 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1038 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1038 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1038 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1038 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1211 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1211 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1211 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1211 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1384 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1384 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1384 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1384 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1557 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1557 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1557 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1557 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=43 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=43 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=43 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=43 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=130 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=130 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=130 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=130 --sim_duration=600 "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=10 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=10 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=10 --sim_duration=600 "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=10 --sim_duration=600 "



#TcpModNewReno
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=0 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=346 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=519 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=692 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1038 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1211 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1384 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1557 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=43 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=130 --sim_duration=600"
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=10 --sim_duration=600"

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=600 --sim_name="ModNewReno-delta-effect" "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=1200 --sim_name="ModNewReno-BW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=1200 --sim_name="ModNewReno-BW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=1200 --sim_name="ModNewReno-BW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=1200 --sim_name="ModNewReno-BW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1730 --sim_duration=1200 --sim_name="ModNewReno-BW-effect" "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=87 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=173 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "

./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpModNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "
./waf --run "examples/tcp/my_dumbbell_sim6 --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=600 --sim_name="ModNewReno-StochBW-effect" "


./waf --run "examples/tcp/my_matplotlib --flows=5 --TcpType="TcpBbrV2" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=300 --sim_name="computational_inv" "
./waf --run "examples/tcp/my_matplotlib --flows=5 --TcpType="TcpVegas" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=300 --sim_name="computational_inv" "
./waf --run "examples/tcp/my_matplotlib --flows=5 --TcpType="TcpCubic" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=300 --sim_name="computational_inv" "
./waf --run "examples/tcp/my_matplotlib --flows=5 --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=865 --sim_duration=300 --sim_name="computational_inv" "
