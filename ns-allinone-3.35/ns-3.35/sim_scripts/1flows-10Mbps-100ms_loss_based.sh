#!/bin/bash
cd ..
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=700 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpWestwood" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=700 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpYeah" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=700 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpVeno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=700 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpIllinois" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=700 --sim_duration=2000 "

./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1200 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpWestwood" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1200 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpYeah" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1200 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpVeno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1200 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpIllinois" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1200 --sim_duration=2000 "

./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1800 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpWestwood" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1800 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpYeah" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1800 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpVeno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1800 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpIllinois" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=1800 --sim_duration=2000 "

./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpNewReno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=500 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpWestwood" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=500 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpYeah" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=500 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpVeno" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=500 --sim_duration=2000 "
./waf --run "examples/tcp/my-dumbbell --flows=1 --sim_name="AggregPerf" --TcpType="TcpIllinois" --btlBW="10Mbps" --btlDelay="100ms" --queueDiscSize=500 --sim_duration=2000 "
