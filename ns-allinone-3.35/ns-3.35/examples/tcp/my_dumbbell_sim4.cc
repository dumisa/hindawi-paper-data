
#include "ns3/my-dumbbell.h"
#include <deque>
#include <eigen3/Eigen/Eigen>
#include <math.h>
#include <iostream>
#include <iomanip>

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
std::vector<SequenceNumber32>  Last_Source_UnackSeq;
std::vector<double>  Last_Source_UnackSeq_time;
std::vector<double>  Last_Time_InFlight;
std::vector<double>  measured_btlBW;
std::vector<uint64_t>  seg_so_far;

std::vector<std::deque<double>> occupancy;
std::vector<std::deque<double>> arrival_rate;
std::vector<std::deque<double>> ntime;
std::vector<std::deque<double>> CoVs;

std::vector<double> rttValue;
std::vector<double> minRtt;
std::vector<uint32_t> cWndValue;
std::vector<uint32_t> inFlightValue;
//uint32_t AckPkts;


Ptr<OutputStreamWrapper> cWndStream;
Ptr<OutputStreamWrapper> RttStream;
Ptr<OutputStreamWrapper> InFlightStream;
Ptr<OutputStreamWrapper> PktsSentStream; 
Ptr<OutputStreamWrapper> seqNumStream; 

/* Function to sort an array using insertion sort*/
double MedianFilter(std::deque<double> arr)
{
    int i, key, j;
    int n = arr. size ();
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;
 
        /* Move elements of arr[0..i-1], that are
        greater than key, to one position ahead
        of their current position */
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }

    return arr[std::ceil(n/2)];
}


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

//numerical differentiation
std::deque<double> TikhonovNumDiff (std::deque<double> x, std::deque<double> y, int k) 
{
  int n = x.size ();
  if (n<1) n=1;
  double a = x.front (), b = x.back ();
  double dt = (b - a)/double(n);
  //std::cout << " " << std::endl << a << " " << b << " " << n << " " << dt << std::endl;
  double t [n];
  MatrixXd A(n, n);
  MatrixXd I(n, n);
  MatrixXd D1(n, n);
  MatrixXd D2(n,n);
  MatrixXd D(n,n);
  Eigen::VectorXd y_hat(n);
  Eigen::VectorXd y_(n);
  Eigen::VectorXd u_(n);
  MatrixXd LHS(n,n);
  Eigen::VectorXd RHS(n);
  double reg_par = 0.01;

  for (int j=0; j < n; j++)
  {
    t[j] = a + j*dt;
  }

  I.setIdentity ();

  A.setZero ();
  for (int i=0; i < n; i++) 
  {
    for (int j=1; j < n+1; j++) 
    {
      if (x[i] <= t[j-1]) A(i,j-1) = 0;
      else if (t[j-1] < x[i] && x[i] < t[j]) A(i,j-1) = x[i] - t[j-1];
      else if (t[j] <= x[i]) A(i,j-1) = dt;
    }
  } 

  D1.setZero (); 
  D1.diagonal(1).setConstant (1/(2*dt));
  D1.diagonal(-1).setConstant (-1/(2*dt));
  D1(0,0) = -1/(2*dt);
  D1(n-1,n-1) = 1/(2*dt);

  D2.setZero (); 
  D2.diagonal(1).setConstant (1/(dt*dt));
  D2.diagonal(0).setConstant (-2/(dt*dt));
  D2.diagonal(-1).setConstant (1/(dt*dt));

  if (k>2) k=2;
  if (k<0) k=0;
  if (k==2) D = I*D1.transpose ()*D2.transpose ();
  else if (k==1) D = I*D1.transpose ();
  else D = I;
  for (int i=0; i < n; i++) 
  {
    y_hat(i) = y[i] - y[0];
    y_(i) = y[i];
  }   

  LHS = A.transpose ()*A + reg_par*D.transpose ()*D;
  RHS = A.transpose ()*y_hat;

  u_ = LHS.inverse () * RHS;

  std::deque<double> u;
   for (int i=0; i < n; i++) 
   {
     //if ( isfinite(u_(i)) ) u.push_back (u_(i));
     u.push_back (u_(i));
   }   

  std::deque<double> out;
  //out.push_back(MedianFilter(u));
  out.push_back(u[u.size () - 1]);
  double Eu;
  Eu = (A*u_ - y_hat).norm () + reg_par*(D*u_).norm ();
  out.push_back(Eu);

  return out;
}

// Trace congestion window
static void CwndTracer (uint16_t source, uint32_t oldval, uint32_t newval)
{
              *cWndStream->GetStream () << Simulator::Now ().GetSeconds () 
                                     << " S" << source
                                     << " " << newval
                                     << " " << (newval*8/rttValue[source])/1e6
                                     << std::endl;
              cWndValue[source] = newval;

}

void InFlightTracer (uint16_t source, uint32_t old, uint32_t inFlight)
{
  NS_UNUSED (old);
  inFlightValue[source] = inFlight;

  double ar_, occ_; 
  double bdp = minRtt[source]*measured_btlBW[source]/pktsize; //pkts
  double now_time = Simulator::Now ().GetSeconds ();

  if (now_time - Last_Time_InFlight[source] >= rttValue[source])
  {
        ar_ = double(inFlight/pktsize/rttValue[source])/bdp;
        occ_ = double(inFlight/pktsize)/bdp;

        if (occupancy[source].size () < 10)
        {  
          occupancy[source].push_back(double(occ_));
          arrival_rate[source].push_back(ar_);
          ntime[source].push_back (now_time);
        } else 
        {
          occupancy[source].pop_front();
          arrival_rate[source].pop_front();
          ntime[source].pop_front ();

          occupancy[source].push_back(double(occ_));
          arrival_rate[source].push_back(ar_);
          ntime[source].push_back (now_time);
        }

          //std::deque<double> u;
          //u = TikhonovNumDiff (ntime, occupancy, 2);
          double d_occupancy = TikhonovNumDiff (ntime[source], occupancy[source], 2)[0];
          //u = TikhonovNumDiff (ntime, occupancy, 2);
          double d_arrival_rate = TikhonovNumDiff (ntime[source], arrival_rate[source], 2)[0];


          double L = occ_;
          L=L; 

          std::deque<double> inter_arrival_time;
          for (uint8_t i; i < arrival_rate.size (); i++)
          {
            inter_arrival_time.push_back (1/arrival_rate[source][i]);
          }

          double CoV = coefficientOfVariation(inter_arrival_time, inter_arrival_time.size ());

          if (CoVs[source].size () < 3) {
            if (isfinite(CoV))  CoVs[source].push_back(CoV);
            else if (CoVs[source].size () < 1) CoVs[source].push_back(0);
            else CoVs[source].push_back(CoVs[source][CoVs[source].size () - 1]);
          }
          else {
            CoVs[source].pop_front();
            if (isfinite(CoV))  CoVs[source].push_back(CoV);
            else CoVs[source].push_back(CoVs[source][CoVs[source].size ()-1]);
          }

          CoV = CoVs[source][2];
          double ratio=d_occupancy/d_arrival_rate;
          //double E_ratio=d_occupancy[1]/d_arrival_rate[1];
          //if(isfinite(ratio) && ratio >= 0)
          //if(not(d_arrival_rate == 0) && atio >= 0)
            if(not(isnan(ratio)))
            {
              *InFlightStream->GetStream () << Simulator::Now ().GetSeconds () 
                                    << " S" << source
                                    << std::fixed << std::setprecision(3)
                                    << " " << bdp
                                    << " " << measured_btlBW[source]*8/1e6
                                    << " " << inFlight/pktsize/bdp
                                    << " " << bdp*ar_*pktsize*8/1e6
                                    << " " << pktsize
                                    << " " << minRtt[source]
                                    << " " << rttValue[source]
                                    << " " << ratio/minRtt[source]
                                    << " " << CoV << "--" << CoVs.size ()
                                    << " " << L
                                    << " " << (1.0+L)*(1.0+L)
                                    << std::endl;
              Last_Time_InFlight[source] = Simulator::Now ().GetSeconds ();
            } 
            //else std::cout << "Not a number or less than zero " << std::endl;
  }

  

}

static void RttTracer (uint16_t source, Time oldval, Time newval)
{
  NS_UNUSED (oldval);
  //static bool first = true;
  // *RttStream->GetStream () << Simulator::Now ().GetSeconds () 
  //                          << " S" << source
  //                          << " " << newval.GetSeconds () << std::endl;
  rttValue[source] = newval.GetSeconds ();
/*  if (first) 
    {
      first=false;
      rttValue[source] = newval.GetSeconds ();
      minRtt = rttValue[source];
    }
  else if (minRtt > rttValue[source]) minRtt = rttValue[source];
  if (minRtt < 0.2) minRtt = 0.2;*/
}

void UnackSequenceTracer (uint16_t source, SequenceNumber32 oldValue, SequenceNumber32 newValue)
{
  //std::cout << "UnackSequenceTracer "  << std::endl;

  double now_time = Simulator::Now ().GetSeconds ();
  UnackSeq = newValue;

  if (now_time - Last_Source_UnackSeq_time[source] > 0)
    {      
      // measured_btlBW[source] = (newValue.GetValue () - Last_Source_UnackSeq[source].GetValue ())*pktsize*8/
      //                                       (now_time - Last_Source_UnackSeq_time[source]); //pkts per seconds
      measured_btlBW[source] = (newValue.GetValue () - Last_Source_UnackSeq[source].GetValue ())/
                                      (now_time - Last_Source_UnackSeq_time[source]); //bytes per time

      seg_so_far[source] += (newValue.GetValue () - Last_Source_UnackSeq[source].GetValue ());

      // *seqNumStream->GetStream () << now_time 
      //                           << " S" << source
      //                           << std::fixed << std::setprecision(3)
      //                           << " " << measured_btlBW[source]/1e6
      //                           << " " << minRtt
      //                           << " " << rttValue
      //                           << std::endl;

            Last_Source_UnackSeq_time[source] = now_time;
            Last_Source_UnackSeq[source] = newValue;
    }
}

static void TcpTxDataTracer (uint16_t source, const Ptr<const Packet> packet, 
                                    const TcpHeader& header, const Ptr<const TcpSocketBase> socket)
{
  if (packet->GetSize () > 100) //to make sure we have the correct trigger
    {
      //*PktsSentStream->GetStream () << Simulator::Now ().GetSeconds ();
      pktsize = packet->GetSize ();

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
      minRtt[source] = socket->GetSocketState ()->m_minRtt.GetSeconds ();
    }   

}

void TraceData (dumbbell& dumbbellSim)
{
    AsciiTraceHelper ascii;
    cWndStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Cwnd.dat");
    //seqNumStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Seq.dat");
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
                << std::endl;

  //std::cout << "PerformanceCalculations " << std::endl;
  Simulator::Schedule (Seconds (time_incr), &PerformanceCalculations, PerformanceDataStream, dumbbellSim);

}

int main (int argc, char *argv[])
{
    std::string TcpType = "TcpBbr";
    //double error_p = 0.0;
    std::string btlBW = "10Mbps";
    std::string accessBW = "10Gbps";
    std::string btlDelay = "100ms";
    std::string accessDelay = "0.01ms";
    std::string queueDisc = "FifoQueueDisc";
    uint32_t queueDiscSize = 0; //packets
    double sim_duration = 1000.0;
    std::string sim_name="sim4_tikhonov";
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
    Last_Source_UnackSeq.reserve(flows);
    Last_Source_UnackSeq_time.reserve(flows);
    Last_Time_InFlight.reserve(flows);
    measured_btlBW.reserve(flows);
    seg_so_far.reserve(flows);
    rttValue.reserve(flows);
    minRtt.reserve(flows);
    cWndValue.reserve(flows);
    inFlightValue.reserve(flows);

    std::deque<double> source_int_;
    std::deque<double> source_int_time;
    source_int_.push_back(0.0);
    source_int_time.push_back (Simulator::Now ().GetSeconds ());

    for (uint8_t i = 0; i < flows; ++i)
    {
      AccumTxPkts[i] = 0;
      AccumRetransPkts[i] = 0;
      AccumLostPkts[i] = 0;
      seg_so_far[i] = 0;
      Last_Source_UnackSeq_time[i] = Simulator::Now ().GetSeconds ();
      Last_Time_InFlight[i] = Simulator::Now ().GetSeconds ();
      occupancy.push_back(source_int_);
      arrival_rate.push_back(source_int_);
      ntime.push_back(source_int_time);
      CoVs.push_back(source_int_);
    }

    //std::cout << "Press ENTER to continue " << Last_Source_UnackSeq.size () <<" " << rttValue[0] << std::endl;
    //getchar();

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