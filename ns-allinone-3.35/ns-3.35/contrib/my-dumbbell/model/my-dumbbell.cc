/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */

#include "my-dumbbell.h"

namespace ns3 {

//function for dumbbell class

dumbbell::~dumbbell () {}
dumbbell::dumbbell () {}

dumbbell::dumbbell(uint8_t flows, std::string TcpType,
    std::string btlBW, std::string btlDelay, 
    std::string accessBW, std::string accessDelay, 
    std::string queueDisc, uint32_t queueDiscSize,
    double error_p, double sim_start, double sim_stop)
{

  // Configure the error model
  // Here we use RateErrorModel with packet error rate
  // Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable> ();
  // uv->SetStream (50);
  // RateErrorModel error_model;
  // error_model.SetRandomVariable (uv);
  // error_model.SetUnit (RateErrorModel::ERROR_UNIT_PACKET);
  // error_model.SetRate (error_p);

    queueDisc = "ns3::" + queueDisc;
  // Configure links
    ns3::PointToPointHelper bottleNeckLink;
    bottleNeckLink.SetDeviceAttribute("DataRate",
            ns3::StringValue(btlBW));
    bottleNeckLink.SetChannelAttribute("Delay",
            ns3::StringValue(btlDelay));
    //bottleNeckLink.SetQueue("ns3::DropTailQueue", "MaxSize", StringValue ("1p"));
    // bottleNeckLink.SetDeviceAttribute("ReceiveErrorModel", PointerValue (&error_model));

    ns3::PointToPointHelper accessLink; //left and right links
    accessLink.SetDeviceAttribute("DataRate",
            ns3::StringValue(accessBW));
    accessLink.SetChannelAttribute("Delay",
            ns3::StringValue(accessDelay));

    //This configuration is important to ensure certainity or minimal requirement for uncontrolled queues
    Config::SetDefault ("ns3::DropTailQueue<Packet>::MaxSize", StringValue ("1p"));

    Config::SetDefault ("ns3::TcpL4Protocol::SocketType", StringValue ("ns3::" + TcpType));
    Config::SetDefault ("ns3::TcpSocket::SndBufSize", UintegerValue (6291456));
    Config::SetDefault ("ns3::TcpSocket::RcvBufSize", UintegerValue (6291456));
    Config::SetDefault ("ns3::TcpSocket::InitialCwnd", UintegerValue (10));
    uint32_t delAckCount = 2;
    Config::SetDefault ("ns3::TcpSocket::DelAckCount", UintegerValue (delAckCount));
    Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (1448));  

    // Create the nodes
    m_routers.Create (2);
    m_leftLeaf.Create (flows);
    m_rightLeaf.Create (flows);
    std::cout << "Nodes created" << std::endl;   

    // Add the link connecting routers
    m_routerDevices = bottleNeckLink.Install (m_routers);
    std::cout << "Bottleneck link installed - created routerdevices" << std::endl;  

    // Add the left side links
    for (uint32_t i = 0; i < flows; ++i)
    {
        NetDeviceContainer c = accessLink.Install (m_routers.Get (0),
                                                 m_leftLeaf.Get (i));
        m_leftRouterDevices.Add (c.Get (0));
        m_leftLeafDevices.Add (c.Get (1));
    }
    std::cout << "Senders connected" << std::endl; 

    // Add the right side links
    for (uint32_t i = 0; i < flows; ++i)
    {
        NetDeviceContainer c = accessLink.Install (m_routers.Get (1),
                                                  m_rightLeaf.Get (i));
        m_rightRouterDevices.Add (c.Get (0));
        m_rightLeafDevices.Add (c.Get (1));
    }   
    std::cout << "Receivers connected" << std::endl; 

    // install internet stack on all nodes
    ns3::InternetStackHelper stack;
    //InstallStack(stack);
    stack.Install(m_routers);
    stack.Install(m_leftLeaf);
    stack.Install(m_rightLeaf);

    // Bottleneck link traffic control configuration
    TrafficControlHelper tchPfifo_btl;
    tchPfifo_btl.SetRootQueueDisc (queueDisc, "MaxSize",
                                 StringValue (std::to_string(queueDiscSize) + "p"));   
    /*TrafficControlHelper tchCoDel_btl;
    tchCoDel.SetRootQueueDisc ("ns3::CoDelQueueDisc");
    Config::SetDefault ("ns3::CoDelQueueDisc::MaxSize", StringValue (std::to_string(queueDiscSize) + "p"));*/

    if (queueDiscSize > 0) {m_btl_queueDiscs = tchPfifo_btl.Install (m_routerDevices);}

    // assign ipv4 addresses
    ns3::Ipv4AddressHelper routerIp =
        ns3::Ipv4AddressHelper("10.3.0.0", "255.255.255.0");
    ns3::Ipv4AddressHelper leftIp =
        ns3::Ipv4AddressHelper("10.1.0.0", "255.255.255.0");
    ns3::Ipv4AddressHelper rightIp =
        ns3::Ipv4AddressHelper("10.2.0.0", "255.255.255.0");
    AssignIpv4Addresses(leftIp, rightIp, routerIp);

    //NS_LOG_INFO ("Initialize Global Routing.");
    Ipv4GlobalRoutingHelper::PopulateRoutingTables ();

    uint16_t port = 50000;
    Address sinkLocalAddress (InetSocketAddress (Ipv4Address::GetAny (), port));
    PacketSinkHelper sinkHelper ("ns3::TcpSocketFactory", sinkLocalAddress);

    m_sender_sockets.reserve (LeftCount ());

    for (uint16_t i = 0; i < LeftCount (); i++)
    {
        AddressValue remoteAddress (InetSocketAddress (m_rightLeafInterfaces.GetAddress (i, 0), port));
        Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (1448));
        BulkSendHelper ftp ("ns3::TcpSocketFactory", Address ());
        ftp.SetAttribute ("Remote", remoteAddress);
        ftp.SetAttribute ("SendSize", UintegerValue (1448));
        ftp.SetAttribute ("MaxBytes", UintegerValue (0));

        ApplicationContainer sourceApp = ftp.Install (m_leftLeaf.Get (i));
        sourceApp.Start (Seconds (sim_start * i));
        sourceApp.Stop (Seconds (sim_stop - 3));

        sinkHelper.SetAttribute ("Protocol", TypeIdValue (TcpSocketFactory::GetTypeId ()));
        ApplicationContainer sinkApp = sinkHelper.Install (m_rightLeaf.Get (i));
        sinkApp.Start (Seconds (sim_start * i));
        sinkApp.Stop (Seconds (sim_stop));
    }
}


void dumbbell::InstallStack (InternetStackHelper stack)
{
  stack.Install (m_routers);
  stack.Install (m_leftLeaf);
  stack.Install (m_rightLeaf);
  std::cout << "Internet stack installed " << std::endl;
}

void dumbbell::AssignIpv4Addresses (Ipv4AddressHelper leftIp,
                  Ipv4AddressHelper rightIp, Ipv4AddressHelper routerIp)
{
  // Assign the router network
  m_routerInterfaces = routerIp.Assign (m_routerDevices);
  // Assign to left side 
  for (uint32_t i = 0; i < LeftCount (); ++i)
    {
      NetDeviceContainer ndc;
      ndc.Add (m_leftLeafDevices.Get (i));
      ndc.Add (m_leftRouterDevices.Get (i));
      Ipv4InterfaceContainer ifc = leftIp.Assign (ndc);
      m_leftLeafInterfaces.Add (ifc.Get (0));
      m_leftRouterInterfaces.Add (ifc.Get (1));
      //std::cout << GetLeftLeafIpv4Address (i) << " " << m_leftRouterInterfaces.GetAddress (i) << std::endl; 
      leftIp.NewNetwork ();
    }
  // Assign to right side
  for (uint32_t i = 0; i < RightCount (); ++i)
    {
      NetDeviceContainer ndc;
      ndc.Add (m_rightLeafDevices.Get (i));
      ndc.Add (m_rightRouterDevices.Get (i));
      Ipv4InterfaceContainer ifc = rightIp.Assign (ndc);
      m_rightLeafInterfaces.Add (ifc.Get (0));
      m_rightRouterInterfaces.Add (ifc.Get (1));
      //std::cout << m_rightLeafInterfaces.GetAddress (i) << " " << m_rightRouterInterfaces.GetAddress (i) << std::endl; 
      rightIp.NewNetwork ();
    }

    std::cout << "IP addresses assigned " << std::endl;
}


void dumbbell::printTopologyConfirmation ()
{
    std::cout << "Confirmation: " << std::endl
              << "     " << RoutersCount () << " routers, " << std::endl
              << "     " << LeftCount ()  << " senders and " << std::endl
              << "     " << RightCount ()  << " receivers" << std::endl
              << "     " << RouterDevicesCount ()  
                        << " BtlNeck Router Devices (netcard)" << std::endl 
              << "     " << LeftRouterDevicesCount ()  
                        << " Left Router Devices (netcard)" << std::endl 
              << "     " << LeftLeafDevicesCount ()  
                        << " Left Leaf Devices (netcard)" << std::endl 
              << "     " << RightRouterDevicesCount ()  
                        << " Right Router Devices (netcard)" << std::endl 
              << "     " << RightLeafDevicesCount ()  << " Right Leaf Devices (netcard)" << std::endl; 

    std::cout << "Done confirmation" << std::endl;
}

void dumbbell::printAssignConfirmation ()
{
    uint32_t flows = LeftCount ();
    std::cout << "Routers IP Addresses:" << std::endl;

    for (uint32_t i=0; i < 2; i++)
    {
        std::cout << "   " << GetRouterIpv4Address (i) << std::endl;
    }

    std::cout << "Left Leaf IP Addresses:" << std::endl;
    for (uint32_t i=0; i < flows; i++)
    {
        std::cout << "   " << GetLeftLeafIpv4Address (i) << std::endl;
    }    

    std::cout << "Left Router IP Addresses:" << std::endl;
    for (uint32_t i=0; i<flows; i++)
    {
        std::cout << "   " << GetLeftRouterIpv4Address (i) << std::endl;
    }

    std::cout << "Right Leaf IP Addresses:" << std::endl;
    for (uint32_t i=0; i<flows; i++)
    {
        std::cout << "   " << GetRightLeafIpv4Address (i) << std::endl;
    }

    std::cout << "Right Router IP Addresses:" << std::endl;
    for (uint32_t i=0; i<flows; i++)
    {
        std::cout << "   " << GetRightRouterIpv4Address (i) << std::endl;
    }
    std::cout << "Done confirmation" << std::endl;
}


}

