
#include "ns3/my-dumbbell.h"

using namespace ns3;

std::string dir;
std::string file_prefix;
double AggregAverageCwndInPkts = 0;
double AggregAveragePktsInFlight = 0;
double AggregAverageRtt =0 ;
uint64_t AccumAggregTxPkts = 0;
uint64_t AccumAggregRetransPkts = 0;
uint64_t AccumAggregRetransPkts2 = 0;
uint64_t AccumAggregLostPkts = 0;
uint64_t AccumAggregLostPkts2 = 0;
uint32_t pktsize = 1448; //1448 default but may be changed during simulation
SequenceNumber32 UnackSeq;

// Trace congestion window
static void CwndTracer (uint32_t oldval, uint32_t newval)
{
  AggregAverageCwndInPkts = 0.5*AggregAverageCwndInPkts + 0.5*newval/double(pktsize); 
}

void TraceCwnd (dumbbell& dumbbellSim)
{
    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace cwnd for Source " << source << " or Node " << nodeId << std::endl;
        //AsciiTraceHelper ascii;
        //Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/" + file_prefix + "cwnd-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/CongestionWindow", MakeCallback (&CwndTracer));
    }
}

void InFlightTracer (uint32_t old, uint32_t inFlight)
{
  NS_UNUSED (old);
  //*TraceInFlightStream->GetStream () << Simulator::Now ().GetSeconds () << " " << inFlight/1448 << std::endl;
  AggregAveragePktsInFlight = 0.5*AggregAveragePktsInFlight + 0.5*inFlight/1448; 
}

void TraceInFlight (dumbbell& dumbbellSim)
{
    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace InFlight Pkts for Source " << source << " or Node " << nodeId << std::endl;
        
        //AsciiTraceHelper ascii;
        //Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/PktsInFlight-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/BytesInFlight", MakeCallback (&InFlightTracer));
    }
}


static void RttTracer (Time oldval, Time newval)
{
  NS_UNUSED (oldval);
  //*rttStream->GetStream () << Simulator::Now ().GetSeconds () << " " << newval.GetSeconds () << std::endl;
  AggregAverageRtt = 0.5*AggregAverageRtt + 0.5*newval.GetSeconds (); 
}

void TraceRtt (dumbbell& dumbbellSim)
{
    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace RTT for Source " << source << " or Node " << nodeId << std::endl;
        
        //AsciiTraceHelper ascii;
        //Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/RTT-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/RTT", MakeCallback (&RttTracer));
    }
}

void UnackSequenceTracer (SequenceNumber32 oldValue, SequenceNumber32 newValue)
{
  //std::cout << "UnackSequenceTracer "  << std::endl;
  UnackSeq = newValue;
}

static void TcpTxDataTracer (const Ptr<const Packet> packet, 
                  const TcpHeader& header, const Ptr<const TcpSocketBase> socket)
{
  //std::cout << "TcpTxDataTracer "  << std::endl;

  if (packet->GetSize () > 100) //to make sure we have the correct trigger
    {
      AccumAggregTxPkts += 1;
      if (socket->GetTxBuffer ()->IsHeadRetransmitted ()) 
        AccumAggregRetransPkts += 1;
      if (socket->GetTxBuffer ()->IsLost(UnackSeq)) 
        AccumAggregLostPkts += 1;
      pktsize = packet->GetSize ();
      AggregAveragePktsInFlight = 0.5*AggregAveragePktsInFlight 
                              + 0.5*socket->GetTxBuffer ()->BytesInFlight ()/double(pktsize);
    }

}

void TraceTcpTxBuffer (dumbbell& dumbbellSim)
{

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();

        std::cout << "Trace TcpSocketBase and unack seq for Source " << source << " or Node " << nodeId 
                  << std::endl;
        
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/TxBuffer/UnackSequence", MakeCallback (&UnackSequenceTracer));
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/Tx", MakeCallback (&TcpTxDataTracer));

    }

}

void TrackProgress(double duration)
{
  static uint32_t count = 0;
  double incr = 10;
  if (count == 0) {std::cout << "Simulation progress " << std::endl;}
  uint32_t perc_duration = 100*Simulator::Now ().GetSeconds ()/duration;
  std::cout << perc_duration << "%" << " " << "....." << std::flush;
  count += 1;
  if (count%5 == 0) {std::cout << std::endl;}

  Simulator::Schedule (Seconds (incr), &TrackProgress, duration);
}

void PerformanceCalculations(Ptr<OutputStreamWrapper> PerformanceDataStream, dumbbell dumbbellSim)
{

  double time_incr;
  //int flows = dumbbellSim.LeftCount ();
  time_incr = 2;

  static bool FirstTime = true;
  if (FirstTime)
  {
        time_incr = time_incr - Simulator::Now ().GetSeconds ();
        FirstTime = false;
        *PerformanceDataStream->GetStream () << "0-Time" << " "
                    << "1-AggregAverageCwndInPkts " << "2-AggregAveragePktsInFlight "
                    << "3-AggregAverageRtt "
                    << "4-AccumAgregTxPkts " << "5-AccumAggregRetransPkts " << "6-AccumAggregLostPkts "
                    << "7-AggregAverageThroughput "
                    << "8-AggregAveragePlr ";
                    //<< "PktSize"; //to confirm if the pkt size is correct
          *PerformanceDataStream->GetStream () << std::endl;
          std::cout << "PerformanceCalculations headings done" << std::endl;
  }

 double AverageGoodput = 0; 
 double AveragePlr = 0;

 if (Simulator::Now ().GetSeconds () > 0)
      AverageGoodput = (1e-6)*(AccumAggregTxPkts - AccumAggregRetransPkts)*8*1448/
                           Simulator::Now ().GetSeconds ();
 if (AccumAggregTxPkts > 0)
      AveragePlr = double(AccumAggregLostPkts)/double(AccumAggregTxPkts);

 if (floor(AveragePlr*1000)/1000 > 0)
    AveragePlr = floor(AveragePlr*10000)/10000;
  else if (floor(AveragePlr*10000)/10000 > 0)
    AveragePlr = floor(AveragePlr*1000000)/1000000;

  *PerformanceDataStream->GetStream () << Simulator::Now ().GetSeconds ()  
                << " " << floor(AggregAverageCwndInPkts) << " " << floor(AggregAveragePktsInFlight)
                << " " << floor(AggregAverageRtt*10000)/10000
                << " " << AccumAggregTxPkts << " " << AccumAggregRetransPkts << " " << AccumAggregLostPkts
                << " " << floor(AverageGoodput*100)/100 
                << " " << AveragePlr;
                //<< " " << pktsize; //to confirm if the pkt size is correct

  *PerformanceDataStream->GetStream () << std::endl;

  //std::cout << "PerformanceCalculations " << std::endl;
  Simulator::Schedule (Seconds (time_incr), &PerformanceCalculations, PerformanceDataStream, dumbbellSim);

}

int main (int argc, char *argv[])
{
    std::string TcpType = "TcpBbr";
    //double error_p = 0.0;
    std::string btlBW = "10Mbps";
    std::string accessBW = "1Gbps";
    std::string btlDelay = "100ms";
    std::string accessDelay = "0.01ms";
    std::string queueDisc = "FifoQueueDisc";
    uint32_t queueDiscSize = 0; //packets
    double sim_duration = 1000.0;
    std::string sim_name="sim_default";
    //bool pcaptracing = false;
  	uint8_t flows = 3;

    CommandLine cmd (__FILE__);
    cmd.AddValue ("TcpType", "Transport protocol to use: TcpNewReno, TcpLinuxReno, "
                "TcpHybla, TcpHighSpeed, TcpHtcp, TcpVegas, TcpScalable, TcpVeno, "
                "TcpBic, TcpYeah, TcpIllinois, TcpWestwood, TcpWestwoodPlus, TcpLedbat, "
                "TcpLp, TcpDctcp, TcpCubic, TcpBbr, TcpBbrV2", TcpType);
    //cmd.AddValue ("error_p", "Packet error rate", error_p);
    cmd.AddValue ("btlBW", "Bottleneck bandwidth", btlBW);
    cmd.AddValue ("btlDelay", "Bottleneck delay", btlDelay);
    cmd.AddValue ("accessBW", "Access link bandwidth", accessBW);
    cmd.AddValue ("accessDelay", "Access link delay", accessDelay);
    cmd.AddValue ("sim_name", "Prefix of output trace file", sim_name);
    cmd.AddValue ("queueDisc", "Queue Disc Type", queueDisc);
    cmd.AddValue ("queueDiscSize", "Size of the Queue Disc", queueDiscSize);
    cmd.AddValue ("flows", "Number of flows", flows);
    cmd.AddValue ("sim_duration", "simulation duration in seconds", sim_duration);
    //cmd.AddValue ("run", "Run index (for setting repeatable seeds)", run);
    //cmd.AddValue ("flow_monitor", "Enable flow monitor", flow_monitor);
    //cmd.AddValue ("pcap_tracing", "Enable or disable PCAP tracing", pcap);
    cmd.Parse (argc, argv);

/*    AccumTxPkts.reserve(flows);
    AccumRetransPkts.reserve(flows);
    AccumLostPkts.reserve(flows);
    AverageCwndInPkts.reserve(flows);
    AveragePktsInFlight.reserve(flows);
    AverageRtt.reserve(flows);
    UnackSeq.reserve(flows);
    PktsInFlight.reserve(flows);
    Rtt.reserve(flows);
    for (uint8_t i = 0; i < flows; ++i)
    {
      AccumTxPkts[i] = 0;
      AccumRetransPkts[i] = 0;
      AccumLostPkts[i] = 0;
      AverageCwndInPkts[i] = 0.0;
      AveragePktsInFlight[i] = 0.0;
      AverageRtt[i] = 0.0;
      PktsInFlight[i] = 0;
      Rtt[i] = 0;
    }*/

    dumbbell dumbbellSim(flows, TcpType, btlBW, btlDelay, accessBW, accessDelay, 
                queueDisc, queueDiscSize, 0, sim_duration);
    ///dumbbellSim.printTopologyConfirmation ();
    //std::cout << "Press ENTER to continue" << std::endl;
    //getchar();
    //dumbbellSim.printAssignConfirmation ();
    //std::cout << "Press ENTER to continue" << std::endl;
    //getchar();

    Simulator::Schedule (Seconds (0.0001), &TraceCwnd, dumbbellSim);
    //Simulator::Schedule (Seconds (0.0002), &TraceInFlight, dumbbellSim);
    Simulator::Schedule (Seconds (0.0003), &TraceRtt, dumbbellSim);
    Simulator::Schedule (Seconds (0.001), &TraceTcpTxBuffer, dumbbellSim);

    //std::cout << "Btl QueueDisks: " << dumbbellSim.GetBtlQueueDiscCount () << std::endl;
    //std::cout << "Btl queue size : " << dumbbellSim.GetBtlQueueDisc (0)->GetMaxSize () << std::endl;
    //std::cout << "Btl queue size : " << dumbbellSim.GetBtlQueueDisc (1)->GetMaxSize () << std::endl;


    // Simulator::Schedule (Seconds (0.00001), &TraceInFlight, prefix_file_name + "-inflight.data");
    // Create a new directory to store the output of the program
    file_prefix = TcpType + "-" + std::to_string(flows) + "-flows-" 
              + btlBW + "-" + btlDelay + "-" 
              + std::to_string(queueDiscSize) + "p-";
    dir = "results/" + sim_name + "/" + std::to_string(flows) + "-flows/" 
              + btlBW + "-" + btlDelay + "/"
              + std::to_string(queueDiscSize) + "p-btlqueue/";
    std::string dirToSave = "mkdir -p " + dir;
    if (system (dirToSave.c_str ()) == -1)
    {
      exit (1);
    } 

    AsciiTraceHelper ascii_file;
    std::string perfData_file_name = file_prefix + "btlqueue-PerfData.dat";
    Ptr<OutputStreamWrapper> PerformanceDataStream = ascii_file.CreateFileStream (dir + perfData_file_name);
    Simulator::Schedule (Seconds (0.001), &PerformanceCalculations, PerformanceDataStream, dumbbellSim);
    Simulator::Schedule (Seconds (0.01), &TrackProgress, sim_duration);


    Simulator::Stop (Seconds (sim_duration));
    Simulator::Run ();
}