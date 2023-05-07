
#include "ns3/my-dumbbell.h"
#include <deque>
#include <eigen3/Eigen/Eigen>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

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
std::vector<uint32_t> cWndValue;
std::vector<uint32_t> inFlightValue;
std::vector<double> rttValue;

std::deque<double> occupancy;
std::deque<double> aggreg_occupancy;
std::deque<double> arrival_rate;
std::deque<double> aggreg_arrival_rate;
std::deque<double> ntime;
std::deque<double> CoVs;
double AverageRttValue=0;
double Last_Time_InFlight=0;

std::vector<uint64_t>  Last_Source_UnackSeq_value;
std::vector<double>  Last_Source_UnackSeq_time;
std::vector<double>  measured_source_btlBW;
double measured_btlBW;

std::default_random_engine generator;
std::poisson_distribution<int> distribution(10);

//uint32_t AckPkts;


Ptr<OutputStreamWrapper> cWndStream;
Ptr<OutputStreamWrapper> RttStream;
Ptr<OutputStreamWrapper> InFlightStream;
Ptr<OutputStreamWrapper> PktsSentStream; 
Ptr<OutputStreamWrapper> seqNumStream; 
Ptr<OutputStreamWrapper> PacketsInQueueStream;

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
void CwndTracer (uint16_t source, uint32_t oldval, uint32_t newval)
{
  //std::cout << "Start CwndTracer " << std::endl;
              // *cWndStream->GetStream () << Simulator::Now ().GetSeconds () 
              //                        << " S" << source
              //                        << " " << newval
              //                        << " " << (newval*8/rttValue[source])/1e6
              //                        << std::endl;
              cWndValue[source] = newval;
  //std::cout << "Complete CwndTracer " << std::endl;

}

void InFlightTracer (uint16_t source, uint32_t old, uint32_t inFlight)
{
  NS_UNUSED (old);
  uint32_t aggregInFlight=0;

  //std::cout << "Start InFlightTracer " << std::endl;

  static bool FirstTime = true;
  if (FirstTime)
  {
        FirstTime = false;
        *InFlightStream->GetStream () << "0-Time"
                    << " 1-Source"
                    << " 2-bdp_pkts"
                    << " 3-arr_Mbps "
                    << " 4-aggreg_ar_Mbps "
                    << " 5-AverageRttValue"
                    << " 6-pktsize"
                    << " 7-mupsi"
                    << " 8-aggreg_mupsi"
                    << " 9-occ_bdp_per_source"
                    << " 10-mupsi_mm1_per_source"
                    << " 11-aggreg_occ_bdp"
                    << " 12-aggreg_mupsi_mm1"
                    << " 13-measured_btlBW_Mbps"
                    << std::endl;

          std::cout << "Inflight headings done" << std::endl;
  }

  inFlightValue[source] = inFlight;

  for (uint32_t i = 0; i < inFlightValue.size (); i++)
  {
    aggregInFlight += inFlightValue[i];
    //std::cout << "inFlight:" << aggregInFlight << std::endl;
  }


  double ar_, occ_, aggreg_ar_, aggreg_occ_; 
  //double bdp = (0.2*10e6/8)/pktsize; //pkts
  double bdp = (0.2*measured_btlBW/8)/pktsize; //pkts
  double now_time = Simulator::Now ().GetSeconds ();

  if (now_time - Last_Time_InFlight >= 0.2)
  {
        ar_ = ((double(inFlight)/double(pktsize)))/AverageRttValue;
        occ_ = (double(inFlight)/double(pktsize));
        aggreg_ar_ = ((double(aggregInFlight)/double(pktsize)))/AverageRttValue;
        aggreg_occ_ = (double(aggregInFlight)/double(pktsize));
        //std::cout << "inFlight:" << aggregInFlight << " ar_: " << AverageRttValue << " occ_:" << occ_ << std::endl;

        if (occupancy.size () < 10)
        {  
          occupancy.push_back(double(occ_));
          arrival_rate.push_back(ar_);
          aggreg_occupancy.push_back(double(aggreg_occ_));
          aggreg_arrival_rate.push_back(aggreg_ar_);
          ntime.push_back (now_time);
        } else 
        {
          occupancy.pop_front();
          arrival_rate.pop_front();
          aggreg_occupancy.pop_front();
          aggreg_arrival_rate.pop_front();
          ntime.pop_front ();

          occupancy.push_back(double(occ_));
          arrival_rate.push_back(ar_);
          aggreg_occupancy.push_back(double(occ_));
          aggreg_arrival_rate.push_back(ar_);
          ntime.push_back (now_time);
        }

          //std::deque<double> u;
          //u = TikhonovNumDiff (ntime, occupancy, 2);
          double d_occupancy = TikhonovNumDiff (ntime, occupancy, 2)[0];
          double aggreg_d_occupancy = TikhonovNumDiff (ntime, aggreg_occupancy, 2)[0];
          //u = TikhonovNumDiff (ntime, occupancy, 2);
          double d_arrival_rate = TikhonovNumDiff (ntime, arrival_rate, 2)[0];
          double aggreg_d_arrival_rate = TikhonovNumDiff (ntime, aggreg_arrival_rate, 2)[0];

          std::deque<double> inter_arrival_time;
          for (uint8_t i=0; i < arrival_rate.size (); i++)
          {
            inter_arrival_time.push_back (1/arrival_rate[i]);
          }

          double CoV = coefficientOfVariation(inter_arrival_time, inter_arrival_time.size ());

          if (CoVs.size () < 3) {
            if (isfinite(CoV))  CoVs.push_back(CoV);
            else if (CoVs.size () < 1) CoVs.push_back(0);
            else CoVs.push_back(CoVs[CoVs.size () - 1]);
          }
          else {
            CoVs.pop_front();
            if (isfinite(CoV))  CoVs.push_back(CoV);
            else CoVs.push_back(CoVs[CoVs.size ()-1]);
          }

          CoV = CoVs[2];
          double psi = d_occupancy/d_arrival_rate;
          double aggreg_occ__psi = aggreg_d_occupancy/aggreg_d_arrival_rate;
          aggreg_occ__psi = aggreg_occ__psi;
          //double E_psi=d_occupancy[1]/d_arrival_rate[1];
          //if(isfinite(psi) && psi >= 0)
          //if(not(d_arrival_rate == 0) && atio >= 0)
            if(not(isnan(psi)) && psi > 0)
            {
              *InFlightStream->GetStream () << Simulator::Now ().GetSeconds () 
                                    << " S" << source
                                    << std::fixed << std::setprecision(6)
                                    << " " << bdp
                                    << " " << ar_*pktsize*8*1e-6
                                    << " " << aggreg_ar_*pktsize*8*1e-6
                                    << " " << AverageRttValue
                                    << " " << pktsize
                                    << " " << psi/0.2
                                    << " " << aggreg_occ__psi/0.2
                                    << " " << occ_/bdp
                                    << " " << pow((1.0+occ_/bdp),2)
                                    << " " << aggreg_occ_/bdp
                                    << " " << pow(1+aggreg_occ_/bdp,2)
                                    << " " << measured_btlBW*1e-6
                                    << std::endl;
              Last_Time_InFlight = Simulator::Now ().GetSeconds ();
            } 
            //else std::cout << "Not a number or less than zero " << std::endl;
  }

  inFlightValue[source] = inFlight;

}

void PacketsInQueueTracer (uint16_t router, uint32_t old, uint32_t newval)
{
  //std::cout << "Start PacketsInQueueTracer " << std::endl;
        // *PacketsInQueueStream->GetStream () << Simulator::Now ().GetSeconds () 
        //                           << " R" << router
        //                           << " " << newval
        //                           << std::endl;
}

void RttTracer (uint16_t source, Time oldval, Time newval)
{
  NS_UNUSED (oldval);
  //std::cout << "Start RttTracer " << std::endl;
   // *RttStream->GetStream () << Simulator::Now ().GetSeconds () 
   //                          << " S" << source
   //                          << " " << newval.GetSeconds () << std::endl;
  double sumRtt=0;
  for (uint32_t i=0; i < rttValue.size (); i++)
  {
    sumRtt += rttValue[i];
  }
  AverageRttValue = sumRtt/rttValue.size ();

  rttValue[source] = newval.GetSeconds ();
  //std::cout << "Complete RttTracer " << std::endl;

}

void UnackSequenceTracer (uint16_t source, SequenceNumber32 oldValue, SequenceNumber32 newValue)
{
  //std::cout << "Start UnackSequenceTracer " << std::endl;
  double now_time = Simulator::Now ().GetSeconds ();
  UnackSeq = newValue;
  static uint64_t BytesSentSoFar=0, lastBytesSoFar=0;
  static double last_time = 0;

  BytesSentSoFar += newValue.GetValue () - oldValue.GetValue ();

  if (now_time - last_time > 0)
    {      
      measured_btlBW = double((BytesSentSoFar - lastBytesSoFar)*8)/(now_time - last_time);
      last_time = now_time;
      lastBytesSoFar = BytesSentSoFar;
    }
    //std::cout << "Complete UnackSequenceTracer " << std::endl;

}

static void TcpTxDataTracer (uint16_t source, const Ptr<const Packet> packet, 
                                    const TcpHeader& header, const Ptr<const TcpSocketBase> socket)
{
  //std::cout << "Start DataTracer " << std::endl;
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
    }  

    //std::cout << "Complete DataTracer " << std::endl; 

}

void TraceData (dumbbell& dumbbellSim)
{
    AsciiTraceHelper ascii;
    //cWndStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Cwnd.dat");
    //seqNumStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Seq.dat");
    //RttStream = ascii.CreateFileStream (dir + "/" + file_prefix + "Rtt.dat");
    InFlightStream = ascii.CreateFileStream (dir + "/" + file_prefix + "InFlight.dat");
    // PktsSentStream = ascii.CreateFileStream (dir + "/" + file_prefix + "PktsSent.dat"); 
    //PacketsInQueueStream = ascii.CreateFileStream (dir + "/" + file_prefix + "PacketsInQueue.dat");

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

    uint32_t nodeId = dumbbellSim.GetLeftRouter ()->GetId ();
    std::cout << "Trace PktsInQueue in router " << 0 << " or Node " << nodeId << std::endl;
    Config::ConnectWithoutContext ("/NodeList/" + std::to_string(nodeId) 
            + "/$ns3::TrafficControlLayer/RootQueueDiscList/0/BytesInQueue", MakeBoundCallback (&PacketsInQueueTracer, nodeId));
    
    std::cout << "Complete setting traces " << std::endl;
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
  static double last_time=0; 
  static uint32_t last_TxPkts=0;
  double now_time;
  double inst_throughput;

  now_time = Simulator::Now ().GetSeconds ();

  static bool FirstTime = true;
  if (FirstTime)
  {
        time_incr = time_incr - now_time;
        FirstTime = false;
        *PerformanceDataStream->GetStream () << "0-Time" << " "
                    << " 1-AverageGoodput "
                    << " 2-AverageThroughput "
                    << " 3-AveragePlr "
                    << " 4-AccumAggregTxPkts "
                    << " 5-AccumAggregLostPkts "
                    << " 6-AccumRetransPkts "
                    << " 7-Instantanous Throughput "
                    << std::endl;

          std::cout << "PerformanceCalculations headings done" << std::endl;
  }
 double AverageGoodput = 0, AverageThroughput = 0; 
 double AveragePlr = 0;

 if (now_time > 0)
 {      
      AverageGoodput = (1e-6)*(AccumAggregTxPkts - AccumAggregRetransPkts)*8*1448/
                           now_time;
      AverageThroughput = (1e-6)*(AccumAggregTxPkts)*8*1448/
                           now_time;
  }
 if (AccumAggregTxPkts > 0)
      AveragePlr = double(AccumAggregRetransPkts)/double(AccumAggregTxPkts);

 if (floor(AveragePlr*1000)/1000 > 0)
    AveragePlr = floor(AveragePlr*10000)/10000;
  else if (floor(AveragePlr*10000)/10000 > 0)
    AveragePlr = floor(AveragePlr*1000000)/1000000;

  if (now_time - last_time > 0)
  {
    inst_throughput = (AccumAggregTxPkts - last_TxPkts)*pktsize*8/(now_time - last_time);
    last_TxPkts = AccumAggregTxPkts;
    last_time = now_time;
  }

  *PerformanceDataStream->GetStream () << Simulator::Now ().GetSeconds ()  
                << " " << floor(AverageGoodput*1000)/1000 
                << " " << floor(AverageThroughput*1000)/1000 
                << " " << AveragePlr
                << " " << AccumAggregTxPkts
                << " " << AccumAggregLostPkts
                << " " << AccumAggregRetransPkts
                << " " << inst_throughput*1e-6
                << std::endl;

  //std::cout << "PerformanceCalculations " << std::endl;
  Simulator::Schedule (Seconds (time_incr), &PerformanceCalculations, PerformanceDataStream, dumbbellSim);

}

void ChangeBW (dumbbell dumbbellSim, uint64_t bps)
{
  PointToPointNetDevice *device=static_cast<PointToPointNetDevice *>(PeekPointer(dumbbellSim.GetRouterDevices (0)));
  device->SetDataRate(DataRate(bps));
}

void StochBW (dumbbell dumbbellSim)
{

  ChangeBW(dumbbellSim, distribution(generator)*1e6);
  Simulator::Schedule (Seconds (0.5), &StochBW, dumbbellSim);
}

int main (int argc, char *argv[])
{
    std::string TcpType = "TcpBbr";
    double error_p = 0.01;
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
    cmd.AddValue ("error_p", "Packet error rate", error_p);
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



    for (uint8_t i = 0; i < flows; ++i)
    {
      AccumTxPkts.push_back (0);
      AccumRetransPkts.push_back (0);
      AccumLostPkts.push_back (0);
      rttValue.push_back (0.2);
      inFlightValue.push_back (0);
      cWndValue.push_back (0);
      Last_Source_UnackSeq_time.push_back (0);
      Last_Source_UnackSeq_value.push_back (0);
      measured_source_btlBW.push_back (0);
    }

    //std::cout << "Press ENTER to continue " << " " <<" " << rttValue[0] << std::endl;
    //getchar();

    dumbbell dumbbellSim(flows, TcpType, btlBW, btlDelay, accessBW, accessDelay, 
                queueDisc, queueDiscSize, error_p, 0, sim_duration);
    ///dumbbellSim.printTopologyConfirmation ();
    //std::cout << "Press ENTER to continue" << std::endl;
    //getchar();
    //dumbbellSim.printAssignConfirmation ();


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
              + std::to_string(queueDiscSize) + "p-"; //2.1-2.01-";
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
    Simulator::Schedule (Seconds (200), &ChangeBW, dumbbellSim, 30e6);
    Simulator::Schedule (Seconds (500), &ChangeBW, dumbbellSim, 10e6);
    Simulator::Schedule (Seconds (800), &ChangeBW, dumbbellSim, 30e6);

    Simulator::Stop (Seconds (sim_duration));
    Simulator::Run ();
}