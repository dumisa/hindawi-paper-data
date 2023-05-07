
#include "ns3/my-dumbbell.h"

using namespace ns3;

std::string dir;
std::string file_prefix;
uint64_t  AccumAggregTxPkts = 0;
uint64_t  AccumAggregRetransPkts = 0;
uint64_t  AccumAggregLostPkts = 0;
uint32_t  pktsize = 1448; //1448 default but may be changed during simulation
SequenceNumber32          UnackSeq;

std::vector<uint64_t>          AccumTxPkts;
std::vector<uint64_t>          AccumRetransPkts;
std::vector<uint64_t>          AccumLostPkts;
std::vector<SequenceNumber32>  Source_UnackSeq;

Ptr<OutputStreamWrapper> AggregTcpTxStream;


// Trace congestion window
static void CwndTracer (Ptr<OutputStreamWrapper> cWndStream, uint32_t oldval, uint32_t newval)
{
  *cWndStream->GetStream () << Simulator::Now ().GetSeconds () << " " << newval/1448 << std::endl;
}

void InFlightTracer (Ptr<OutputStreamWrapper> InFlightStream, uint32_t old, uint32_t inFlight)
{
  NS_UNUSED (old);
  *InFlightStream->GetStream () << Simulator::Now ().GetSeconds () << " " << inFlight/1448 << std::endl;
}

static void RttTracer (Ptr<OutputStreamWrapper> RttStream, Time oldval, Time newval)
{
  NS_UNUSED (oldval);
  *RttStream->GetStream () << Simulator::Now ().GetSeconds () << " " << newval.GetSeconds () << std::endl;
}


void UnackSequenceTracer (uint16_t source, SequenceNumber32 oldValue, SequenceNumber32 newValue)
{
  //std::cout << "UnackSequenceTracer "  << std::endl;
  UnackSeq = newValue;
  Source_UnackSeq[source] = newValue;
}

static void TcpTxDataTracer (Ptr<OutputStreamWrapper> SourceTcpTxStream, uint16_t source, const Ptr<const Packet> packet, 
                                    const TcpHeader& header, const Ptr<const TcpSocketBase> socket)
{

  if (packet->GetSize () > 100) //to make sure we have the correct trigger
    {
      *AggregTcpTxStream->GetStream () << Simulator::Now ().GetSeconds ();
      *SourceTcpTxStream->GetStream () << Simulator::Now ().GetSeconds ();

      AccumAggregTxPkts += 1;
      AccumTxPkts[source] += 1;
      *AggregTcpTxStream->GetStream () << " " << 1;
      *SourceTcpTxStream->GetStream () << " " << 1;
      if (socket->GetTxBuffer ()->IsHeadRetransmitted ()) 
      {
        AccumAggregRetransPkts += 1;
        AccumRetransPkts[source] += 1;
        *AggregTcpTxStream->GetStream () << " " << 1;
        *SourceTcpTxStream->GetStream () << " " << 1;
      }
      else
      {
        *AggregTcpTxStream->GetStream () << " " << 0;
        *SourceTcpTxStream->GetStream () << " " << 0;
      }
      if (socket->GetTxBuffer ()->IsLost(UnackSeq)) 
      {
        AccumAggregLostPkts += 1;
        AccumLostPkts[source] += 1;
      }

      *AggregTcpTxStream->GetStream () << std::endl;
      *SourceTcpTxStream->GetStream () << std::endl;
    }

}



void TraceData (dumbbell& dumbbellSim)
{
    AsciiTraceHelper ascii;
    AggregTcpTxStream = ascii.CreateFileStream (dir + "/" + file_prefix + "AggregPktsSent" + ".dat");   

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace cwnd for Source " << source << " or Node " << nodeId << std::endl;

        Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/" + file_prefix + "Cwnd-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/CongestionWindow", MakeBoundCallback (&CwndTracer, stream));
    }


    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace RTT for Source " << source << " or Node " << nodeId << std::endl;

        Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/" + file_prefix + "Rtt-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/RTT", MakeBoundCallback (&RttTracer, stream));
    }

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace InFlight for Source " << source << " or Node " << nodeId << std::endl;

        Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/" + file_prefix + "InFlight-S" + std::to_string(source) + ".dat");
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/BytesInFlight", MakeBoundCallback (&InFlightTracer, stream));
    }

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace TxPkts for Source " << source << " or Node " << nodeId << std::endl;

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/TxBuffer/UnackSequence", MakeBoundCallback (&UnackSequenceTracer, source));

        Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (dir + "/" + file_prefix + "PktsSent-S" + std::to_string(source) + ".dat");   
        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/Tx", MakeBoundCallback (&TcpTxDataTracer, stream, source));

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
                    << "1-AccumAgregTxPkts " << "2-AccumAggregRetransPkts " << "3-AccumAggregLostPkts "
                    << "4-AggregAverageThroughput "
                    << "5-AggregAveragePlr "
                    << std::endl;

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
                << " " << AccumAggregTxPkts << " " << AccumAggregRetransPkts << " " << AccumAggregLostPkts
                << " " << floor(AverageGoodput*1000)/1000 
                << " " << AveragePlr
                << std::endl;

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

    AccumTxPkts.reserve(flows);
    AccumRetransPkts.reserve(flows);
    AccumLostPkts.reserve(flows);
    Source_UnackSeq.reserve(flows);

    for (uint8_t i = 0; i < flows; ++i)
    {
      AccumTxPkts[i] = 0;
      AccumRetransPkts[i] = 0;
      AccumLostPkts[i] = 0;
    }

    dumbbell dumbbellSim(flows, TcpType, btlBW, btlDelay, accessBW, accessDelay, 
                queueDisc, queueDiscSize, 0, sim_duration);
    ///dumbbellSim.printTopologyConfirmation ();
    //std::cout << "Press ENTER to continue" << std::endl;
    //getchar();
    //dumbbellSim.printAssignConfirmation ();
    //std::cout << "Press ENTER to continue" << std::endl;
    //getchar();

    Simulator::Schedule (Seconds (0.0001), &TraceData, dumbbellSim);
    //Simulator::Schedule (Seconds (0.0002), &TraceInFlight, dumbbellSim);
/*    Simulator::Schedule (Seconds (0.0003), &TraceRtt, dumbbellSim);
    Simulator::Schedule (Seconds (0.001), &TraceTcpTxBuffer, dumbbellSim);*/

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