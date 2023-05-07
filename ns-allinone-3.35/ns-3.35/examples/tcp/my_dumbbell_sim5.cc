
#include "ns3/my-dumbbell.h"
#include <deque>
#include <eigen3/Eigen/Eigen>
#include <math.h>

using namespace ns3;
using Eigen::MatrixXd;

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
std::vector<uint64_t>          PktsInFlight;
std::vector<SequenceNumber32>  Source_UnackSeq;

std::deque<double> occupancy;
std::deque<double> arrival_rate;
std::deque<double> ntime;
double rttValue;
double minRtt;
uint64_t inFlightValue;
double measured_btlBW;
uint32_t AckPkts;


Ptr<OutputStreamWrapper> cWndStream;
Ptr<OutputStreamWrapper> RttStream;
Ptr<OutputStreamWrapper> InFlightStream;
Ptr<OutputStreamWrapper> PktsSentStream; 

// Function to find mean of given array.
double mean(std::deque<double> arr, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + arr[i];
    return sum / n;
}

// Function to find standard deviation
// of given array.
double standardDeviation(std::deque<double> arr, int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + (arr[i] - mean(arr, n)) *
                    (arr[i] - mean(arr, n));
 
    return sqrt(sum / (n - 1));
}

// Function to find coefficient of variation.
double coefficientOfVariation(std::deque<double> arr, int n)
{
   return standardDeviation(arr, n) / mean(arr, n);
}

// Trace congestion window
static void CwndTracer (uint16_t source, uint32_t oldval, uint32_t newval)
{
  //static double lastAggTxPkts = 0; 
  static double last_time = Simulator::Now ().GetSeconds ();
  double now_time = Simulator::Now ().GetSeconds ();

  if (now_time - last_time > 2*minRtt)
  {
        //ar_ = 1448*(AccumAggregTxPkts - AccumAggregRetransPkts - last_total_acks)/
                  (Simulator::Now ().GetSeconds () - last_time);
        //ar_ = 1448.0*(AccumAggregTxPkts - AccumAggregRetransPkts - lastAggTxPkts)/(now_time - last_time);

        if (occupancy.size () < 10)
        {  
          occupancy.push_back (double(inFlightValue));
          arrival_rate.push_back(double(newval)/rttValue);
          ntime.push_back (now_time);
        } else 
        {
          occupancy.pop_front ();
          arrival_rate.pop_front ();
          ntime.pop_front ();

          occupancy.push_back (double(inFlightValue));
          arrival_rate.push_back(double(newval)/rttValue);
          ntime.push_back (now_time);
        }

        double CoV;

        CoV = coefficientOfVariation(arrival_rate, arrival_rate.size ());

        double L;
        double bdp = (minRtt*10e6/8);
        L = double(inFlightValue)/bdp; 

        double dt = now_time - last_time;
        double num = (occupancy[occupancy.size ()-1] - occupancy[occupancy.size () - 3])/(2*dt);
        double den = (arrival_rate[arrival_rate.size ()-1] -  arrival_rate[ arrival_rate.size () - 3])/(2*dt);
        //double ratio = (dt/2) * (occupancy[occupancy.size ()-1] - occupancy[occupancy.size () - 3]) / 
          //                (occupancy[occupancy.size ()-1] - 2*occupancy[occupancy.size () - 2] + occupancy[occupancy.size () - 3]);
        double ratio = num/den;

        if ( (ratio >= 1) && not isinf(ratio) && not isnan(ratio))
          {
              *cWndStream->GetStream () << Simulator::Now ().GetSeconds () 
                                    //<< " S" << source
                                    << " " << newval/1448
                                    << " " << measured_btlBW*8/1e6
                                    << " " << rttValue
                                    << " " << ratio/minRtt
                                    << " " << CoV
                                    << " " << L
                                    << " " << bdp/1448.0
                                    << " " << (1.0+L)*(1.0+L) << std::endl;
              last_time = Simulator::Now ().GetSeconds ();
          }
            //else std::cout << "Not a number or less than zero " << std::endl;

        //lastAggTxPkts = AccumAggregTxPkts - AccumAggregRetransPkts;
          //last_time = Simulator::Now ().GetSeconds ();
        

  }

}

void InFlightTracer (uint16_t source, uint32_t old, uint32_t inFlight)
{
  NS_UNUSED (old);
  static double ar_, occ_; 
  //static double lastAggTxPkts=0;
  static double last_time = Simulator::Now ().GetSeconds ();
  double L;
  double bdp = minRtt*10e6/8;

  double now_time = Simulator::Now ().GetSeconds ();
  if (now_time - last_time > 2*minRtt)
  {
        ar_ = double(inFlight/bdp)/rttValue;
        occ_ = double(inFlight/bdp);

        if (occupancy.size () < 10)
        {  
          occupancy.push_back (double(occ_));
          arrival_rate.push_back(ar_);
          ntime.push_back (now_time);
        } else 
        {
          occupancy.pop_front ();
          arrival_rate.pop_front ();
          ntime.pop_front ();

          occupancy.push_back (double(occ_));
          arrival_rate.push_back(ar_);
          ntime.push_back (now_time);
        }

        L = double(inFlight)/bdp; 
        L=L;

          std::deque<double> inter_arrival_time;
          for (uint8_t i; i < occupancy.size (); i++)
          {
            inter_arrival_time.push_back (arrival_rate[i]);
          }

        double CoV = coefficientOfVariation(inter_arrival_time, occupancy.size ());

        double dt = now_time - last_time;
        double num = (occupancy[occupancy.size ()-1] - occupancy[occupancy.size () - 3])/(2*dt);
        double den = (arrival_rate[arrival_rate.size ()-1] -  arrival_rate[ arrival_rate.size () - 3])/(2*dt);
        //double ratio = (dt/2) * (occupancy[occupancy.size ()-1] - occupancy[occupancy.size () - 3]) / 
          //                (occupancy[occupancy.size ()-1] - 2*occupancy[occupancy.size () - 2] + occupancy[occupancy.size () - 3]);
        double ratio = num/den;

        if ( (ratio >= 1) && not isinf(ratio) && not isnan(ratio))
            {
              *InFlightStream->GetStream () << Simulator::Now ().GetSeconds () 
                                    //<< " S" << source
                                    //<< " " << measured_btlBW
                                    //<< " " << occupancy.size ()
                                    << " " << ar_
                                    << " " << rttValue
                                    //<< " " << d_occupancy[1]
                                    //<< " " << d_arrival_rate[1]
                                    << " " << ratio/minRtt
                                    //<< " " << E_ratio/minRtt
                                    << " " << CoV
                                    << " " << L
                                    << " " << (1.0+L)*(1.0+L) << std::endl;
              last_time = Simulator::Now ().GetSeconds ();
            } 
            //else std::cout << "Not a number or less than zero " << std::endl;
    }


  inFlightValue = inFlight;
}

static void RttTracer (uint16_t source, Time oldval, Time newval)
{
  NS_UNUSED (oldval);
  static bool first = true;
  // *RttStream->GetStream () << Simulator::Now ().GetSeconds () 
  //                          << " S" << source
  //                          << " " << newval.GetSeconds () << std::endl;
  rttValue = newval.GetSeconds ();
  if (first) 
    {
      first=false;
      rttValue = newval.GetSeconds ();
      minRtt = rttValue;
    }
  else if (minRtt > rttValue) minRtt = rttValue;
  if (minRtt < 0.2) minRtt = 0.2;
}


void UnackSequenceTracer (uint16_t source, SequenceNumber32 oldValue, SequenceNumber32 newValue)
{
  //std::cout << "UnackSequenceTracer "  << std::endl;
  UnackSeq = newValue;
  Source_UnackSeq[source] = newValue;
}

static void TcpTxDataTracer (uint16_t source, const Ptr<const Packet> packet, 
                                    const TcpHeader& header, const Ptr<const TcpSocketBase> socket)
{

  static double last_time = Simulator::Now ().GetSeconds ();

  if (packet->GetSize () > 100) //to make sure we have the correct trigger
    {
      //*PktsSentStream->GetStream () << Simulator::Now ().GetSeconds ();

      AccumAggregTxPkts += 1;
      AccumTxPkts[source] += 1;
      //*PktsSentStream->GetStream () << " S" << source << " " << 1;
      if (socket->GetTxBuffer ()->IsHeadRetransmitted ()) 
      {
        AccumAggregRetransPkts += 1;
        AccumRetransPkts[source] += 1;
        //*PktsSentStream->GetStream () << " " << 1;
      }
      else
      {
        //*PktsSentStream->GetStream () << " " << 0;
      }
      if (socket->GetTxBuffer ()->IsLost(UnackSeq)) 
      {
        AccumAggregLostPkts += 1;
        AccumLostPkts[source] += 1;
      }
      //*PktsSentStream->GetStream () << std::endl;
      if (socket->GetTxBuffer ()->GetSacked () > 0 && 
                                Simulator::Now ().GetSeconds () - last_time > 0) 
        {
          measured_btlBW = double(socket->GetTxBuffer ()->GetSacked ())/
                    (Simulator::Now ().GetSeconds () - last_time);
          last_time = Simulator::Now ().GetSeconds ();
        } else
      {
        //startSentTime = Simulator::Now ().GetSeconds ();
      }

    }


    

}



void TraceData (dumbbell& dumbbellSim)
{
    AsciiTraceHelper ascii;
    cWndStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Cwnd.dat");
    // RttStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Rtt.dat");
    InFlightStream = ascii.CreateFileStream (dir + "/" + file_prefix + "InFlight.dat");
    // PktsSentStream = ascii.CreateFileStream (dir + "/" + file_prefix + "PktsSent.dat");  

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace cwnd for Source " << source << " or Node " << nodeId << std::endl;

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/CongestionWindow", 
              MakeBoundCallback (&CwndTracer, source));
    }


    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace RTT for Source " << source << " or Node " << nodeId << std::endl;

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/RTT", MakeBoundCallback (&RttTracer, source));
    }

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace InFlight for Source " << source << " or Node " << nodeId << std::endl;

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/BytesInFlight", 
              MakeBoundCallback (&InFlightTracer, source));
    }

    for (uint32_t source = 0; source < dumbbellSim.LeftCount (); ++source)
    {
        uint32_t nodeId = dumbbellSim.GetLeftLeaf (source)->GetId ();
        std::cout << "Trace TxPkts for Source " << source << " or Node " << nodeId << std::endl;

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/TxBuffer/UnackSequence", MakeBoundCallback (&UnackSequenceTracer, source));

        Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
              + "/$ns3::TcpL4Protocol/SocketList/*/Tx", MakeBoundCallback (&TcpTxDataTracer, source));

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
                << " " << AckPkts
                << " " << measured_btlBW
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
    std::string sim_name="sim5_normal";
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