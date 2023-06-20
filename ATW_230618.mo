package ATW
  import Modelica.Media;

  model pipe_
    replaceable package Medium = Media.Water.StandardWaterOnePhase;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    import Modelica.Constants.pi;
    // parameters
    parameter SI.Length pipe_length = 5 "pipe length [m]";
    parameter SI.Length pipe_diameterInner = 0.025 "pipe inner diameter [m]";
    parameter SI.Length pipe_diameterOuter = 0.030 "pipe outer diameter [m]";
    parameter SI.ThermalConductivity pipe_TC = 398 "Cupper [W/(m.K)]";
    parameter SI.SpecificHeatCapacity pipe_c = 386 "[J/(kg.K)]";
    parameter SI.Density pipe_rho = 8960 "[kg/m3]";
    parameter SI.CoefficientOfHeatTransfer air_h = 10 "W/(m2.K)";
    parameter SI.Temperature pipe_initialTemp = 273.15 + 30 "pipe initial temperature [K]";
    // fluid (pipe)
    Modelica.Fluid.Pipes.DynamicPipe pipe(redeclare package Medium = Medium, T_start = pipe_initialTemp, allowFlowReversal = false, diameter = pipe_diameterInner, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, height_ab = 0, isCircular = true, length = pipe_length, m_flow_start = 0.0001, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, modelStructure = Modelica.Fluid.Types.ModelStructure.av_b, momentumDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nNodes = 5, nParallel = 1, roughness = 0, use_HeatTransfer = true, use_T_start = true) annotation(
      Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // fluid - wall
    Modelica.Thermal.HeatTransfer.Components.Convection convection_pipeInner[pipe.nNodes] annotation(
      Placement(visible = true, transformation(origin = {-50, 26}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Modelica.Blocks.Sources.Constant const_GcInner[pipe.nNodes](each k = 10000000) annotation(
      Placement(visible = true, transformation(origin = {-90, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // wall
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor[pipe.nNodes](each G = 2*pi*pipe_TC*(pipe_length/pipe.nNodes)/log((pipe_diameterInner + (pipe_diameterOuter - pipe_diameterInner)/2)/pipe_diameterInner)) annotation(
      Placement(visible = true, transformation(origin = {-22, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1[pipe.nNodes](each G = 2*pi*pipe_TC*(pipe_length/pipe.nNodes)/log(pipe_diameterOuter/(pipe_diameterInner + (pipe_diameterOuter - pipe_diameterInner)/2))) annotation(
      Placement(visible = true, transformation(origin = {30, 46}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor[pipe.nNodes](each C = pipe_rho*pi*(pipe_diameterOuter^2 - pipe_diameterInner^2)/4*pipe_length/pipe.nNodes*pipe_c, each T(start = pipe_initialTemp)) annotation(
      Placement(visible = true, transformation(origin = {4, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // wall - air
    Modelica.Thermal.HeatTransfer.Components.Convection convection_pipeOuter[pipe.nNodes] annotation(
      Placement(visible = true, transformation(origin = {66, 46}, extent = {{10, 10}, {-10, -10}}, rotation = -180)));
    Modelica.Blocks.Sources.Constant const_GcOuter[pipe.nNodes](each k = air_h*pi*pipe_diameterOuter*pipe_length/pipe.nNodes) annotation(
      Placement(visible = true, transformation(origin = {44, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // air
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature[pipe.nNodes](each T = 273.15 - 10) annotation(
      Placement(visible = true, transformation(origin = {104, 46}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    // ports
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.port_b, port_b) annotation(
      Line(points = {{-40, 0}, {120, 0}}, color = {0, 127, 255}));
    connect(pipe.port_a, port_a) annotation(
      Line(points = {{-60, 0}, {-120, 0}}, color = {0, 127, 255}));
    connect(convection_pipeInner.fluid, pipe.heatPorts) annotation(
      Line(points = {{-50, 16}, {-50, 4}}, color = {191, 0, 0}, thickness = 0.5));
    connect(convection_pipeInner.Gc, const_GcInner.y) annotation(
      Line(points = {{-60, 26}, {-78, 26}}, color = {0, 0, 127}, thickness = 0.5));
    connect(convection_pipeInner.solid, thermalConductor.port_a) annotation(
      Line(points = {{-50, 36}, {-50, 46}, {-32, 46}}, color = {191, 0, 0}, thickness = 0.5));
    connect(thermalConductor.port_b, heatCapacitor.port) annotation(
      Line(points = {{-12, 46}, {4, 46}, {4, 62}}, color = {191, 0, 0}, thickness = 0.5));
    connect(heatCapacitor.port, thermalConductor1.port_b) annotation(
      Line(points = {{4, 62}, {4, 46}, {20, 46}}, color = {191, 0, 0}, thickness = 0.5));
    connect(thermalConductor1.port_a, convection_pipeOuter.solid) annotation(
      Line(points = {{40, 46}, {56, 46}}, color = {191, 0, 0}, thickness = 0.5));
    connect(const_GcOuter.y, convection_pipeOuter.Gc) annotation(
      Line(points = {{56, 76}, {66, 76}, {66, 56}}, color = {0, 0, 127}, thickness = 0.5));
    connect(convection_pipeOuter.fluid, fixedTemperature.port) annotation(
      Line(points = {{76, 46}, {94, 46}}, color = {191, 0, 0}, thickness = 0.5));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 40}, {100, -40}}), Line(origin = {-0.28, -70.7614}, points = {{-59.7236, 5.20272}, {60.2764, 5.20272}, {40.2764, 15.2027}, {60.2764, 5.20272}, {40.2764, -4.79728}}, color = {67, 134, 234}, thickness = 0.5), Rectangle(origin = {0, 45}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 5}, {100, -5}}), Rectangle(origin = {0, -45}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 5}, {100, -5}}), Line(origin = {3.31729, 79.4407}, points = {{-13.3308, -28.9993}, {8.66919, -4.99929}, {-15.3308, -4.99929}, {14.6692, 29.0007}}, color = {170, 0, 0}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 8), Text(origin = {-3, -2}, extent = {{-85, 24}, {85, -24}}, textString = "%name", textStyle = {TextStyle.Bold})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-140, 100}, {140, -20}})));
  end pipe_;

  model test_pipe_
    replaceable package Medium = Media.Water.StandardWaterOnePhase;
    ATW.pipe_ pipe_(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, T = 273.15 + 80, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-54, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 5, height = 0.167, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-86, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature_a(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-22, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature_b(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {24, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe_.port_b, boundary.ports[1]) annotation(
      Line(points = {{10, 0}, {40, 0}}, color = {0, 127, 255}));
    connect(pipe_.port_a, boundary1.ports[1]) annotation(
      Line(points = {{-10, 0}, {-36, 0}}, color = {0, 127, 255}));
    connect(ramp.y, boundary1.m_flow_in) annotation(
      Line(points = {{-75, 8}, {-56, 8}}, color = {0, 0, 127}));
    connect(pipe_.port_a, temperature_a.port) annotation(
      Line(points = {{-10, 0}, {-22, 0}, {-22, 8}}, color = {0, 127, 255}));
    connect(pipe_.port_b, temperature_b.port) annotation(
      Line(points = {{10, 0}, {24, 0}, {24, 8}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3));
  end test_pipe_;

  model HEX_HP
    replaceable package Medium = Media.Water.StandardWaterOnePhase;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-128, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter Real K = 150 "";
    parameter Real T = 500 "";
    parameter SI.Temperature InitialTemp = 273.15 + 30 "pipe initial temperature [K]";
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = InitialTemp, V = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, use_HeatTransfer = true, use_T_start = true, use_portsData = false, nPorts = 2) annotation(
      Placement(visible = true, transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput freq annotation(
      Placement(visible = true, transformation(origin = {-126, 46}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-98, 90}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
      Placement(visible = true, transformation(origin = {-50, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.FirstOrder firstOrder(T = T, initType = Modelica.Blocks.Types.Init.InitialOutput, k = K, y_start = 0) annotation(
      Placement(visible = true, transformation(origin = {-82, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(port_a, volume.ports[1]) annotation(
      Line(points = {{-120, 0}, {0, 0}}));
    connect(port_b, volume.ports[2]) annotation(
      Line(points = {{120, 0}, {0, 0}}));
    connect(prescribedHeatFlow.port, volume.heatPort) annotation(
      Line(points = {{-40, 46}, {-26, 46}, {-26, 10}, {-10, 10}}, color = {191, 0, 0}));
    connect(freq, firstOrder.u) annotation(
      Line(points = {{-126, 46}, {-94, 46}}, color = {0, 0, 127}));
    connect(firstOrder.y, prescribedHeatFlow.Q_flow) annotation(
      Line(points = {{-70, 46}, {-60, 46}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 40}, {100, -40}}), Line(origin = {-0.28, -70.7614}, points = {{-59.7236, 5.20272}, {60.2764, 5.20272}, {40.2764, 15.2027}, {60.2764, 5.20272}, {40.2764, -4.79728}}, color = {67, 134, 234}, thickness = 0.5), Rectangle(origin = {0, 45}, fillColor = {180, 180, 180}, fillPattern = FillPattern.Solid, extent = {{-100, 5}, {100, -5}}), Rectangle(origin = {0, -45}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 5}, {100, -5}}), Text(origin = {-3, -2}, extent = {{-85, 24}, {85, -24}}, textString = "%name", textStyle = {TextStyle.Bold}), Rectangle(origin = {0, 90}, fillColor = {255, 176, 123}, fillPattern = FillPattern.Solid, extent = {{-100, 40}, {100, -40}}), Line(origin = {5.7, 49.8}, points = {{-13.3308, -28.9993}, {8.66919, -4.99929}, {-15.3308, -4.99929}, {14.6692, 29.0007}}, color = {170, 0, 0}, thickness = 1, arrow = {Arrow.Open, Arrow.Open}, arrowSize = 8), Rectangle(origin = {0, 135}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 5}, {100, -5}})}, coordinateSystem(extent = {{-120, 140}, {120, -80}})),
      Diagram(coordinateSystem(extent = {{-140, 100}, {140, -20}})));
  end HEX_HP;

  model test_HEX_HP
    replaceable package Medium = Media.Water.StandardWaterOnePhase;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-54, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {50, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, T = 273.15 + 80, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 0.167, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-86, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature_a(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-22, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature_b(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {24, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.HEX_HP hex_hp(K = 2, redeclare package Medium = Medium, T = 100) annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-12, 14}, {12, 36}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = 100, offset = 0, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-32, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ramp.y, boundary1.m_flow_in) annotation(
      Line(points = {{-75, 8}, {-56, 8}}, color = {0, 0, 127}));
    connect(boundary1.ports[1], hex_hp.port_a) annotation(
      Line(points = {{-36, 0}, {-10, 0}}, color = {0, 127, 255}));
    connect(boundary.ports[1], hex_hp.port_b) annotation(
      Line(points = {{40, 0}, {10, 0}}, color = {0, 127, 255}));
    connect(hex_hp.port_a, temperature_a.port) annotation(
      Line(points = {{-10, 0}, {-22, 0}, {-22, 8}}, color = {0, 127, 255}));
    connect(hex_hp.port_b, temperature_b.port) annotation(
      Line(points = {{10, 0}, {24, 0}, {24, 8}}, color = {0, 127, 255}));
    connect(ramp1.y, hex_hp.freq) annotation(
      Line(points = {{-20, 68}, {-14, 68}, {-14, 10}, {-10, 10}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
  end test_HEX_HP;

  model room
    //  replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    replaceable package Medium = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-90, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Length room_W = 3.6 "Width [m]";
    parameter SI.Length room_D = 7.2 "Depth [m]";
    parameter SI.Length room_H = 0.4 "Height [m]";
    parameter SI.Length wall_thick = 0.1 "wall thicks [m]";
    parameter SI.Area room_area = 2*(room_W*room_D + room_D*room_H + room_H*room_W);
    parameter SI.Volume room_volume = room_W*room_D*room_H;
    parameter SI.Volume wall_volume = room_area*wall_thick;
    parameter SI.Temperature room_initialTemp = 273.15 + 20 "room initial temperature [K]";
    parameter SI.Temperature wall_initialTemp = 273.15 + 15 "room initial temperature [K]";
    parameter SI.Density wall_density = 144 "[kg/m3]";
    parameter SI.SpecificHeatCapacity wall_specHeatCap = 1168 "[J/kg.K]";
    parameter SI.ThermalConductivity wall_TC = 0.3 "[W/(m.K)]";
    parameter SI.CoefficientOfHeatTransfer air_h = 10 "???[W/(m2.K)]";
    parameter SI.Temperature T_amb = 273.15 - 5 "room initial temperature [K]";
    parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
    // def
    Modelica.Fluid.Vessels.ClosedVolume air(redeclare package Medium = Medium, T_start = room_initialTemp, V = room_volume, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_HeatTransfer = true, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation(
      Placement(visible = true, transformation(origin = {-38, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection_room annotation(
      Placement(visible = true, transformation(origin = {-2, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = wall_volume*wall_density*wall_specHeatCap, T(fixed = true, start = wall_initialTemp)) annotation(
      Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection_amb annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T = T_amb) annotation(
      Placement(visible = true, transformation(origin = {130, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const_room(k = air_h*room_area) annotation(
      Placement(visible = true, transformation(origin = {-26, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = wall_TC*room_area*2/wall_thick) annotation(
      Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = wall_TC*room_area*2/wall_thick) annotation(
      Placement(visible = true, transformation(origin = {74, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = Medium, redeclare package Medium2 = Medium, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 500, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
      Placement(visible = true, transformation(origin = {-24, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pumpRA(redeclare package Medium = Medium, N_nominal = 1000, T_start = T_amb, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true, use_powerCharacteristic = false) annotation(
      Placement(visible = true, transformation(origin = {38, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 10, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {16, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pumpSA(redeclare package Medium = Medium, N_nominal = 1000, T_start = T_amb, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true, use_powerCharacteristic = false) annotation(
      Placement(visible = true, transformation(origin = {72, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    ATW.pipe_ pipeRA(redeclare package Medium = Medium, pipe_initialTemp = T_amb) annotation(
      Placement(visible = true, transformation(origin = {6, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.pipe_ pipeSA(redeclare package Medium = Medium, pipe_initialTemp = T_amb, pipe_length = 0.01) annotation(
      Placement(visible = true, transformation(origin = {6, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT amb(redeclare package Medium = Medium, T = T_amb, nPorts = 2, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {134, -68}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperatureSA(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-76, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //
    // Interfaces
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_rad annotation(
      Placement(visible = true, transformation(origin = {-100, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_conv annotation(
      Placement(visible = true, transformation(origin = {-100, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput T annotation(
      Placement(visible = true, transformation(origin = {146, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperatureRA(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-78, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperatureOA(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {38, -114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperatureOA1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {82, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(port_rad, air.heatPort) annotation(
      Line(points = {{-100, 28}, {-62, 28}, {-62, 0}, {-52, 0}}, color = {191, 0, 0}));
    connect(port_conv, air.heatPort) annotation(
      Line(points = {{-100, -18}, {-62, -18}, {-62, 0}, {-52, 0}}, color = {191, 0, 0}));
    connect(air.heatPort, temperatureSensor.port) annotation(
      Line(points = {{-52, 0}, {-56, 0}, {-56, -94}, {-48, -94}}, color = {191, 0, 0}));
    connect(air.heatPort, convection_room.fluid) annotation(
      Line(points = {{-52, 0}, {-54, 0}, {-54, -14}, {-20, -14}, {-20, 0}, {-12, 0}}, color = {191, 0, 0}));
    connect(convection_amb.fluid, fixedTemperature.port) annotation(
      Line(points = {{110, 0}, {122, 0}}, color = {191, 0, 0}));
    connect(const_room.y, convection_room.Gc) annotation(
      Line(points = {{-15, 52}, {-2, 52}, {-2, 10}}, color = {0, 0, 127}));
    connect(convection_room.solid, thermalConductor_wall1.port_a) annotation(
      Line(points = {{8, 0}, {16, 0}}, color = {191, 0, 0}));
    connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
      Line(points = {{36, 0}, {50, 0}, {50, 20}}, color = {191, 0, 0}));
    connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
      Line(points = {{50, 20}, {50, 0}, {64, 0}}, color = {191, 0, 0}));
    connect(thermalConductor_wall2.port_a, convection_amb.solid) annotation(
      Line(points = {{84, 0}, {90, 0}}, color = {191, 0, 0}));
    connect(const_room.y, convection_amb.Gc) annotation(
      Line(points = {{-14, 52}, {100, 52}, {100, 10}}, color = {0, 0, 127}));
    connect(temperatureSensor.T, T) annotation(
      Line(points = {{-27, -94}, {146, -94}}, color = {0, 0, 127}));
    connect(air.ports[1], Lossnay.port_a1) annotation(
      Line(points = {{-42, -10}, {-36, -10}, {-36, -60}, {-34, -60}}, color = {0, 127, 255}));
    connect(Lossnay.port_b2, air.ports[2]) annotation(
      Line(points = {{-34, -72}, {-40, -72}, {-40, -10}, {-42, -10}}, color = {0, 127, 255}));
    connect(ramp.y, pumpRA.N_in) annotation(
      Line(points = {{28, -28}, {39, -28}, {39, -50}, {38, -50}}, color = {0, 0, 127}));
    connect(ramp.y, pumpSA.N_in) annotation(
      Line(points = {{28, -28}, {72, -28}, {72, -62}}, color = {0, 0, 127}));
    connect(Lossnay.port_b1, pipeRA.port_a) annotation(
      Line(points = {{-14, -60}, {-4, -60}}, color = {0, 127, 255}));
    connect(pipeRA.port_b, pumpRA.port_a) annotation(
      Line(points = {{16, -60}, {28, -60}}, color = {0, 127, 255}));
    connect(Lossnay.port_a2, pipeSA.port_b) annotation(
      Line(points = {{-14, -72}, {-4, -72}}, color = {0, 127, 255}));
    connect(pipeSA.port_a, pumpSA.port_b) annotation(
      Line(points = {{16, -72}, {62, -72}}, color = {0, 127, 255}));
    connect(pumpRA.port_b, amb.ports[1]) annotation(
      Line(points = {{48, -60}, {124, -60}, {124, -68}}, color = {0, 127, 255}));
    connect(pumpSA.port_a, amb.ports[2]) annotation(
      Line(points = {{82, -72}, {124, -72}, {124, -68}}, color = {0, 127, 255}));
    connect(Lossnay.port_b2, temperatureSA.port) annotation(
      Line(points = {{-34, -72}, {-55, -72}, {-55, -70}, {-76, -70}}, color = {0, 127, 255}));
    connect(temperatureRA.port, Lossnay.port_a1) annotation(
      Line(points = {{-78, -36}, {-55, -36}, {-55, -60}, {-34, -60}}, color = {0, 127, 255}));
    connect(Lossnay.port_a2, temperatureOA.port) annotation(
      Line(points = {{-14, -72}, {-14, -124}, {38, -124}}, color = {0, 127, 255}));
    connect(pipeSA.port_a, temperatureOA1.port) annotation(
      Line(points = {{16, -72}, {50, -72}, {50, -118}, {82, -118}}, color = {0, 127, 255}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 170, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -3}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Line(origin = {19.9746, 97.1102}, points = {{-62.0003, -38.0003}, {5.99968, 25.9997}, {-4.00032, -28.0003}, {61.9997, 37.9997}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Open}, arrowSize = 15)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-120, 80}, {160, -100}})));
  end room;

  model test_room
    replaceable package Medium = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1000, offset = 0, startTime = 5) annotation(
      Placement(visible = true, transformation(origin = {-62, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
      Placement(visible = true, transformation(origin = {-34, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.room room(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {34, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(prescribedHeatFlow.port, room.port_conv) annotation(
      Line(points = {{-24, 30}, {0, 30}, {0, 32}, {24, 32}}, color = {191, 0, 0}));
    connect(ramp.y, prescribedHeatFlow.Q_flow) annotation(
      Line(points = {{-50, 30}, {-44, 30}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
  end test_room;

  model roomCycle
    replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-116, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Temperature Water_initT = 273.15 + 30;
    parameter SI.Temperature Rad_initT = 273.15 + 40;
    ATW.pipe_ pipe_1(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-102, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.pipe_ pipe_2(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {6, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.pipe_ pipe_3(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {6, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    ATW.pipe_ pipe_4(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-100, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 110, offset = 0, startTime = 5) annotation(
      Placement(visible = true, transformation(origin = {-66, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Medium, N_nominal = 1750, T_start = Water_initT, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.034, 0.04}, head_nominal = {39.0, 27.0, 22.8}), m_flow_start = 0.0001, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, p_a_start = 101325, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower(V_flow_nominal = {0, 0.034, 0.04}, W_nominal = {5000, 14700, 17000}), use_N_in = true, use_T_start = true, use_powerCharacteristic = true) annotation(
      Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = Water_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325) annotation(
      Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    // Interfaces
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Con annotation(
      Placement(visible = true, transformation(origin = {68, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Rad annotation(
      Placement(visible = true, transformation(origin = {68, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(port_a, pipe_1.port_a) annotation(
      Line(points = {{-130, 38}, {-112, 38}}));
    connect(pipe_1.port_b, pump.port_a) annotation(
      Line(points = {{-92, 38}, {-42, 38}}, color = {0, 127, 255}));
    connect(pump.port_b, pipe_2.port_a) annotation(
      Line(points = {{-22, 38}, {-4, 38}}, color = {0, 127, 255}));
    connect(pipe_2.port_b, rad.port_a) annotation(
      Line(points = {{16, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
    connect(rad.port_b, pipe_3.port_a) annotation(
      Line(points = {{24, -4}, {24, -32}, {16, -32}}, color = {0, 127, 255}));
    connect(pipe_3.port_b, pipe_4.port_a) annotation(
      Line(points = {{-4, -32}, {-90, -32}}, color = {0, 127, 255}));
    connect(pipe_4.port_b, port_b) annotation(
      Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
    connect(ramp.y, pump.N_in) annotation(
      Line(points = {{-54, 68}, {-32, 68}, {-32, 48}}, color = {0, 0, 127}));
    connect(rad.heatPortCon, port_Con) annotation(
      Line(points = {{32, 8}, {46, 8}, {46, 22}, {68, 22}}, color = {191, 0, 0}));
    connect(rad.heatPortRad, port_Rad) annotation(
      Line(points = {{32, 4}, {46, 4}, {46, -12}, {68, -12}}, color = {191, 0, 0}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-140, 100}, {80, -40}})));
  end roomCycle;

  model test_roomcycle
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.roomCycle roomCycle(redeclare package Medium = Medium, Water_initT = 273.15 + 30) annotation(
      Placement(visible = true, transformation(origin = {-28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 273.15 + 80, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-94, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-96, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.room room(redeclare package Medium = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {40, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.ports[1], roomCycle.port_a) annotation(
      Line(points = {{-84, 36}, {-38, 36}, {-38, 32}}, color = {0, 127, 255}));
    connect(roomCycle.port_b, boundary1.ports[1]) annotation(
      Line(points = {{-38, 22}, {-76, 22}, {-76, -4}, {-86, -4}}, color = {0, 127, 255}));
    connect(roomCycle.port_Con, room.port_conv) annotation(
      Line(points = {{-18, 32}, {22, 32}, {22, 42}, {30, 42}}, color = {191, 0, 0}));
    connect(roomCycle.port_Rad, room.port_rad) annotation(
      Line(points = {{-18, 22}, {28, 22}, {28, 32}, {30, 32}}, color = {191, 0, 0}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
  end test_roomcycle;

  model plant
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Media.Air.DryAirNasa;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-170, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Volume BT_vol = 0.001 "Buffer Tank Volume [m3]";
    parameter SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
    parameter SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
    parameter SI.Temperature Water_initT = 273.15 + 30;
    //
    ATW.HEX_HP hex_hp(InitialTemp = Water_initT, redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-56, 62}, extent = {{-12, 14}, {12, 36}}, rotation = 90)));
    ATW.pipe_ pipe_4(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-42, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    ATW.pipe_ pipe_1(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-56, 88}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    ATW.pipe_ pipe_2(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {10, 104}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Modelica.Fluid.Vessels.ClosedVolume BT(redeclare package Medium = Medium, T_start = BT_initT, V = BT_vol, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 4, p_start = BT_initP, use_HeatTransfer = false, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {50, 76}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    ATW.pipe_ pipe_3(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {22, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Medium, N_nominal = 1000, T_start = Water_initT, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true) annotation(
      Placement(visible = true, transformation(origin = {-12, -12}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 38.5, offset = 0, startTime = 5) annotation(
      Placement(visible = true, transformation(origin = {-38, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_ph boundary_beg(redeclare package Medium = Medium, nPorts = 1, p = 101325, use_h_in = true) annotation(
      Placement(visible = true, transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Sources.Boundary_pT boundary_end(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-90, 34}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    ATW.room room(redeclare package Medium = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.roomCycle roomCycle(redeclare package Medium = Medium, Rad_initT = Water_initT, Water_initT = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // sensors
    Modelica.Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-90, 6}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-26, 104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // interfaces
    Modelica.Blocks.Interfaces.RealInput f annotation(
      Placement(visible = true, transformation(origin = {-146, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-98, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput T annotation(
      Placement(visible = true, transformation(origin = {154, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(hex_hp.port_b, pipe_1.port_a) annotation(
      Line(points = {{-56, 72}, {-56, 78}}, color = {0, 127, 255}));
    connect(pipe_2.port_b, BT.ports[1]) annotation(
      Line(points = {{20, 104}, {42, 104}, {42, 76}, {40, 76}}, color = {0, 127, 255}));
    connect(BT.ports[2], pipe_3.port_a) annotation(
      Line(points = {{40, 76}, {40, -12}, {32, -12}}, color = {0, 127, 255}));
    connect(pipe_4.port_a, pump.port_b) annotation(
      Line(points = {{-32, -12}, {-22, -12}}, color = {0, 127, 255}));
    connect(pump.port_a, pipe_3.port_b) annotation(
      Line(points = {{-2, -12}, {12, -12}}, color = {0, 127, 255}));
    connect(ramp.y, pump.N_in) annotation(
      Line(points = {{-27, -38}, {-13, -38}, {-13, -22}}, color = {0, 0, 127}));
    connect(boundary_beg.ports[1], hex_hp.port_a) annotation(
      Line(points = {{-56, 44}, {-56, 52}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port_b, boundary_end.ports[1]) annotation(
      Line(points = {{-90, 16}, {-90, 24}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port_a, pipe_4.port_b) annotation(
      Line(points = {{-90, -4}, {-90, -12}, {-52, -12}}, color = {0, 127, 255}));
    connect(specificEnthalpy.h_out, boundary_beg.h_in) annotation(
      Line(points = {{-78, 6}, {-60, 6}, {-60, 22}}, color = {0, 0, 127}));
    connect(f, hex_hp.freq) annotation(
      Line(points = {{-146, 62}, {-76, 62}, {-76, 48}, {-64, 48}, {-64, 52}}, color = {0, 0, 127}));
    connect(pipe_1.port_b, massFlowRate.port_a) annotation(
      Line(points = {{-56, 98}, {-56, 104}, {-36, 104}}, color = {0, 127, 255}));
    connect(massFlowRate.port_b, pipe_2.port_a) annotation(
      Line(points = {{-16, 104}, {0, 104}}, color = {0, 127, 255}));
    connect(roomCycle.port_Con, room.port_conv) annotation(
      Line(points = {{80, 56}, {100, 56}}, color = {191, 0, 0}));
    connect(roomCycle.port_Rad, room.port_rad) annotation(
      Line(points = {{80, 46}, {100, 46}}, color = {191, 0, 0}));
    connect(BT.ports[3], roomCycle.port_b) annotation(
      Line(points = {{40, 76}, {24, 76}, {24, 46}, {60, 46}}, color = {0, 127, 255}));
    connect(BT.ports[4], roomCycle.port_a) annotation(
      Line(points = {{40, 76}, {30, 76}, {30, 56}, {60, 56}}, color = {0, 127, 255}));
    connect(room.T, T) annotation(
      Line(points = {{120, 50}, {154, 50}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-180, 120}, {160, -60}})),
      Icon(graphics = {Rectangle(fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {1, 2}, extent = {{-63, 38}, {63, -38}}, textString = "%name", textStyle = {TextStyle.Bold})}));
  end plant;

  model controller
    extends Modelica.Icons.Example;
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.plant plant(redeclare package Medium = Medium, redeclare package MediumAir = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {-2, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 273.15 + 20) annotation(
      Placement(visible = true, transformation(origin = {-84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback annotation(
      Placement(visible = true, transformation(origin = {-56, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.PI pi(T = 300, initType = Modelica.Blocks.Types.Init.InitialState, k = 5) annotation(
      Placement(visible = true, transformation(origin = {-28, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(const.y, feedback.u1) annotation(
      Line(points = {{-72, 18}, {-64, 18}}, color = {0, 0, 127}));
    connect(plant.T, feedback.u2) annotation(
      Line(points = {{8, 18}, {20, 18}, {20, 0}, {-56, 0}, {-56, 10}}, color = {0, 0, 127}));
    connect(feedback.y, pi.u) annotation(
      Line(points = {{-46, 18}, {-40, 18}}, color = {0, 0, 127}));
    connect(pi.y, plant.f) annotation(
      Line(points = {{-16, 18}, {-12, 18}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-100, 80}, {20, -20}})));
  end controller;

  model room_test
    //  replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    replaceable package Medium = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-86, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Length room_W = 3.6 "Width [m]";
    parameter SI.Length room_D = 7.2 "Depth [m]";
    parameter SI.Length room_H = 0.4 "Height [m]";
    parameter SI.Length wall_thick = 0.1 "wall thicks [m]";
    parameter SI.Area room_area = 2*(room_W*room_D + room_D*room_H + room_H*room_W);
    parameter SI.Volume room_volume = room_W*room_D*room_H;
    parameter SI.Volume wall_volume = room_area*wall_thick;
    parameter SI.Temperature room_initialTemp = 273.15 + 20 "room initial temperature [K]";
    parameter SI.Temperature wall_initialTemp = 273.15 + 15 "room initial temperature [K]";
    parameter SI.Density wall_density = 144 "[kg/m3]";
    parameter SI.SpecificHeatCapacity wall_specHeatCap = 1168 "[J/kg.K]";
    parameter SI.ThermalConductivity wall_TC = 0.3 "[W/(m.K)]";
    parameter SI.CoefficientOfHeatTransfer air_h = 10 "???[W/(m2.K)]";
    parameter SI.Temperature T_amb = 273.15 - 5 "room initial temperature [K]";
    parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
    // def
    Modelica.Fluid.Vessels.ClosedVolume air(redeclare package Medium = Medium, T_start = room_initialTemp, V = room_volume, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = 101325, use_HeatTransfer = true, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = Medium, redeclare package Medium2 = Medium, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 500, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
      Placement(visible = true, transformation(origin = {-24, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //
    // Interfaces
    Modelica.Fluid.Sources.Boundary_pT amb(redeclare package Medium = Medium, T = T_amb, nPorts = 2, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {108, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    pipe_ pipeSA(redeclare package Medium = Medium, pipe_initialTemp = T_amb, pipe_length = 1) annotation(
      Placement(visible = true, transformation(origin = {6, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    pipe_ pipeRA(redeclare package Medium = Medium, pipe_initialTemp = T_amb, pipe_length = 1) annotation(
      Placement(visible = true, transformation(origin = {6, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pumpSA(redeclare package Medium = Medium, N_nominal = 1750, T_start = T_amb, V = 1, allowFlowReversal = false, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.034, 0.04}, head_nominal = {39.0, 27.0, 22.8}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower(V_flow_nominal = {0, 0.034, 0.04}, W_nominal = {5000, 14700, 17000}), use_HeatTransfer = false, use_N_in = true, use_T_start = true, use_powerCharacteristic = true) annotation(
      Placement(visible = true, transformation(origin = {72, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 100, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {16, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pumpRA(redeclare package Medium = Medium, N_nominal = 1750, T_start = T_amb, V = 1, allowFlowReversal = false, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.034, 0.04}, head_nominal = {39.0, 27.0, 22.8}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower(V_flow_nominal = {0, 0.034, 0.04}, W_nominal = {5000, 14700, 17000}), use_HeatTransfer = false, use_N_in = true, use_T_start = true, use_powerCharacteristic = true) annotation(
      Placement(visible = true, transformation(origin = {38, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Lossnay.port_b2, air.ports[1]) annotation(
      Line(points = {{-34, -72}, {-42, -72}, {-42, -10}}, color = {0, 127, 255}));
    connect(Lossnay.port_a2, pipeSA.port_b) annotation(
      Line(points = {{-14, -72}, {-4, -72}}, color = {0, 127, 255}));
    connect(Lossnay.port_b1, pipeRA.port_a) annotation(
      Line(points = {{-14, -60}, {-4, -60}}, color = {0, 127, 255}));
    connect(pumpSA.port_a, amb.ports[1]) annotation(
      Line(points = {{82, -72}, {90, -72}, {90, -70}, {98, -70}}, color = {0, 127, 255}));
    connect(pipeSA.port_a, pumpSA.port_b) annotation(
      Line(points = {{16, -72}, {62, -72}}, color = {0, 127, 255}));
    connect(ramp.y, pumpSA.N_in) annotation(
      Line(points = {{28, -28}, {72, -28}, {72, -62}}, color = {0, 0, 127}));
    connect(pumpRA.port_b, amb.ports[2]) annotation(
      Line(points = {{48, -60}, {90, -60}, {90, -70}, {98, -70}}, color = {0, 127, 255}));
    connect(pipeRA.port_b, pumpRA.port_a) annotation(
      Line(points = {{16, -60}, {28, -60}}, color = {0, 127, 255}));
    connect(ramp.y, pumpRA.N_in) annotation(
      Line(points = {{28, -28}, {39, -28}, {39, -50}, {38, -50}}, color = {0, 0, 127}));
    connect(air.ports[2], Lossnay.port_a1) annotation(
      Line(points = {{-42, -10}, {-40, -10}, {-40, -60}, {-34, -60}}, color = {0, 127, 255}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 170, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -3}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Line(origin = {19.9746, 97.1102}, points = {{-62.0003, -38.0003}, {5.99968, 25.9997}, {-4.00032, -28.0003}, {61.9997, 37.9997}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Open}, arrowSize = 15)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-100, 20}, {120, -100}})),
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
  end room_test;

  model room_test1
    //  replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    replaceable package Medium = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-86, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Length room_W = 3.6 "Width [m]";
    parameter SI.Length room_D = 7.2 "Depth [m]";
    parameter SI.Length room_H = 0.4 "Height [m]";
    parameter SI.Length wall_thick = 0.1 "wall thicks [m]";
    parameter SI.Area room_area = 2*(room_W*room_D + room_D*room_H + room_H*room_W);
    parameter SI.Volume room_volume = room_W*room_D*room_H;
    parameter SI.Volume wall_volume = room_area*wall_thick;
    parameter SI.Temperature room_initialTemp = 273.15 + 20 "room initial temperature [K]";
    parameter SI.Temperature wall_initialTemp = 273.15 + 15 "room initial temperature [K]";
    parameter SI.Density wall_density = 144 "[kg/m3]";
    parameter SI.SpecificHeatCapacity wall_specHeatCap = 1168 "[J/kg.K]";
    parameter SI.ThermalConductivity wall_TC = 0.3 "[W/(m.K)]";
    parameter SI.CoefficientOfHeatTransfer air_h = 10 "???[W/(m2.K)]";
    parameter SI.Temperature T_amb = 273.15 - 5 "room initial temperature [K]";
    parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
    // def
    Modelica.Fluid.Vessels.ClosedVolume air(redeclare package Medium = Medium, T_start = room_initialTemp, V = room_volume, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = 101325, use_HeatTransfer = true, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = Medium, redeclare package Medium2 = Medium, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 500, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
      Placement(visible = true, transformation(origin = {-24, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //
    // Interfaces
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = T_amb, nPorts = 1, p = 101325 - 100) annotation(
      Placement(visible = true, transformation(origin = {12, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = T_amb, nPorts = 2, p = 101325 + 100) annotation(
      Placement(visible = true, transformation(origin = {54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(air.ports[1], Lossnay.port_a1) annotation(
      Line(points = {{-42, -10}, {-40, -10}, {-40, -60}, {-34, -60}}, color = {0, 127, 255}));
    connect(Lossnay.port_b2, air.ports[2]) annotation(
      Line(points = {{-34, -72}, {-42, -72}, {-42, -10}}, color = {0, 127, 255}));
    connect(boundary.ports[1], Lossnay.port_b1) annotation(
      Line(points = {{22, -26}, {31, -26}, {31, -60}, {-14, -60}}, color = {0, 127, 255}));
    connect(boundary1.ports[1], Lossnay.port_a2) annotation(
      Line(points = {{64, -40}, {78, -40}, {78, -72}, {-14, -72}}, color = {0, 127, 255}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 170, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -3}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Line(origin = {19.9746, 97.1102}, points = {{-62.0003, -38.0003}, {5.99968, 25.9997}, {-4.00032, -28.0003}, {61.9997, 37.9997}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Open}, arrowSize = 15)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-100, 20}, {120, -100}})),
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
  end room_test1;

  model room_test2
    //replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    replaceable package Medium = Buildings.Media.Air;
    //  replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-86, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Length room_W = 3.6 "Width [m]";
    parameter SI.Length room_D = 7.2 "Depth [m]";
    parameter SI.Length room_H = 0.4 "Height [m]";
    parameter SI.Length wall_thick = 0.1 "wall thicks [m]";
    parameter SI.Area room_area = 2*(room_W*room_D + room_D*room_H + room_H*room_W);
    parameter SI.Volume room_volume = room_W*room_D*room_H;
    parameter SI.Volume wall_volume = room_area*wall_thick;
    parameter SI.Temperature room_initialTemp = 273.15 + 20 "room initial temperature [K]";
    parameter SI.Temperature wall_initialTemp = 273.15 + 15 "room initial temperature [K]";
    parameter SI.Density wall_density = 144 "[kg/m3]";
    parameter SI.SpecificHeatCapacity wall_specHeatCap = 1168 "[J/kg.K]";
    parameter SI.ThermalConductivity wall_TC = 0.3 "[W/(m.K)]";
    parameter SI.CoefficientOfHeatTransfer air_h = 10 "???[W/(m2.K)]";
    parameter SI.Temperature T_amb = 273.15 + 15 "room initial temperature [K]";
    parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
    // def
    Modelica.Fluid.Vessels.ClosedVolume air(redeclare package Medium = Medium, T_start = room_initialTemp, V = room_volume, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = 101325, use_HeatTransfer = true, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = Medium, redeclare package Medium2 = Medium, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 500, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
      Placement(visible = true, transformation(origin = {-24, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //
    // Interfaces
    Modelica.Fluid.Sources.Boundary_pT amb(redeclare package Medium = Medium, T = T_amb, nPorts = 2, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {108, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    pipe_ pipeSA(redeclare package Medium = Medium, pipe_initialTemp = T_amb, pipe_length = 1) annotation(
      Placement(visible = true, transformation(origin = {6, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    pipe_ pipeRA(redeclare package Medium = Medium, pipe_initialTemp = T_amb, pipe_length = 1) annotation(
      Placement(visible = true, transformation(origin = {6, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 0, offset = 100, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {16, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pumpSA(redeclare package Medium = Medium, N_nominal = 1000, T_start = T_amb, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true, use_powerCharacteristic = false) annotation(
      Placement(visible = true, transformation(origin = {72, -72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-70, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Lossnay.port_b2, air.ports[1]) annotation(
      Line(points = {{-34, -72}, {-42, -72}, {-42, -10}}, color = {0, 127, 255}));
    connect(Lossnay.port_a2, pipeSA.port_b) annotation(
      Line(points = {{-14, -72}, {-4, -72}}, color = {0, 127, 255}));
    connect(Lossnay.port_b1, pipeRA.port_a) annotation(
      Line(points = {{-14, -60}, {-4, -60}}, color = {0, 127, 255}));
    connect(air.ports[2], Lossnay.port_a1) annotation(
      Line(points = {{-42, -10}, {-40, -10}, {-40, -60}, {-34, -60}}, color = {0, 127, 255}));
    connect(ramp.y, pumpSA.N_in) annotation(
      Line(points = {{28, -28}, {72, -28}, {72, -62}}, color = {0, 0, 127}));
    connect(pipeSA.port_a, pumpSA.port_b) annotation(
      Line(points = {{16, -72}, {62, -72}}, color = {0, 127, 255}));
    connect(pumpSA.port_a, amb.ports[1]) annotation(
      Line(points = {{82, -72}, {90, -72}, {90, -70}, {98, -70}}, color = {0, 127, 255}));
    connect(pipeRA.port_b, amb.ports[2]) annotation(
      Line(points = {{16, -60}, {98, -60}, {98, -70}}, color = {0, 127, 255}));
    connect(Lossnay.port_b2, temperature.port) annotation(
      Line(points = {{-34, -72}, {-70, -72}, {-70, -44}}, color = {0, 127, 255}));
    annotation(
      Icon(graphics = {Rectangle(fillColor = {170, 170, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -3}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Line(origin = {19.9746, 97.1102}, points = {{-62.0003, -38.0003}, {5.99968, 25.9997}, {-4.00032, -28.0003}, {61.9997, 37.9997}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Open}, arrowSize = 15)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-100, 20}, {120, -100}})),
      experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.06));
  end room_test2;

  model plant2
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Media.Air.DryAirNasa;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-170, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // parameters
    parameter SI.Volume BT_vol = 0.001 "Buffer Tank Volume [m3]";
    parameter SI.Volume BT_heiht = 0.1 "Buffer Tank height [m]";
    parameter SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
    parameter SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
    parameter SI.Temperature Water_initT = 273.15 + 30;
    //
    ATW.HEX_HP hex_hp(InitialTemp = Water_initT, redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-56, 62}, extent = {{-12, 14}, {12, 36}}, rotation = 90)));
    ATW.pipe_ pipe_4(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-42, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    ATW.pipe_ pipe_1(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {-56, 88}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    ATW.pipe_ pipe_2(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {10, 104}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    ATW.pipe_ pipe_3(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {22, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Medium, N_nominal = 1000, T_start = Water_initT, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true) annotation(
      Placement(visible = true, transformation(origin = {-12, -12}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 38.5, offset = 0, startTime = 5) annotation(
      Placement(visible = true, transformation(origin = {-38, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_ph boundary_beg(redeclare package Medium = Medium, nPorts = 1, p = 101325, use_h_in = true) annotation(
      Placement(visible = true, transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Sources.Boundary_pT boundary_end(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-90, 34}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    ATW.room room(redeclare package Medium = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {112, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ATW.roomCycle roomCycle(redeclare package Medium = Medium, Rad_initT = Water_initT, Water_initT = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {72, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // sensors
    Modelica.Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-90, 6}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-26, 104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // interfaces
    Modelica.Blocks.Interfaces.RealInput f annotation(
      Placement(visible = true, transformation(origin = {-146, 62}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-98, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput T annotation(
      Placement(visible = true, transformation(origin = {156, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.Storage.Stratified BT(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_heiht, nSeg = 5, p_start = BT_initP) annotation(
      Placement(visible = true, transformation(origin = {42, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(hex_hp.port_b, pipe_1.port_a) annotation(
      Line(points = {{-56, 72}, {-56, 78}}, color = {0, 127, 255}));
    connect(pipe_4.port_a, pump.port_b) annotation(
      Line(points = {{-32, -12}, {-22, -12}}, color = {0, 127, 255}));
    connect(pump.port_a, pipe_3.port_b) annotation(
      Line(points = {{-2, -12}, {12, -12}}, color = {0, 127, 255}));
    connect(ramp.y, pump.N_in) annotation(
      Line(points = {{-27, -38}, {-13, -38}, {-13, -22}}, color = {0, 0, 127}));
    connect(boundary_beg.ports[1], hex_hp.port_a) annotation(
      Line(points = {{-56, 44}, {-56, 52}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port_b, boundary_end.ports[1]) annotation(
      Line(points = {{-90, 16}, {-90, 24}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port_a, pipe_4.port_b) annotation(
      Line(points = {{-90, -4}, {-90, -12}, {-52, -12}}, color = {0, 127, 255}));
    connect(specificEnthalpy.h_out, boundary_beg.h_in) annotation(
      Line(points = {{-78, 6}, {-60, 6}, {-60, 22}}, color = {0, 0, 127}));
    connect(f, hex_hp.freq) annotation(
      Line(points = {{-146, 62}, {-76, 62}, {-76, 48}, {-64, 48}, {-64, 52}}, color = {0, 0, 127}));
    connect(pipe_1.port_b, massFlowRate.port_a) annotation(
      Line(points = {{-56, 98}, {-56, 104}, {-36, 104}}, color = {0, 127, 255}));
    connect(massFlowRate.port_b, pipe_2.port_a) annotation(
      Line(points = {{-16, 104}, {0, 104}}, color = {0, 127, 255}));
    connect(roomCycle.port_Con, room.port_conv) annotation(
      Line(points = {{82, 91}, {102, 91}}, color = {191, 0, 0}));
    connect(roomCycle.port_Rad, room.port_rad) annotation(
      Line(points = {{82, 81}, {102, 81}}, color = {191, 0, 0}));
    connect(room.T, T) annotation(
      Line(points = {{122.2, 86}, {156.2, 86}}, color = {0, 0, 127}));
  connect(pipe_2.port_b, BT.port_a) annotation(
      Line(points = {{20, 104}, {32, 104}, {32, 70}}, color = {0, 127, 255}));
  connect(BT.port_b, pipe_3.port_a) annotation(
      Line(points = {{52, 70}, {56, 70}, {56, -12}, {32, -12}}, color = {0, 127, 255}));
  connect(BT.fluPorVol[1], roomCycle.port_a) annotation(
      Line(points = {{40, 70}, {40, 92}, {62, 92}}, color = {0, 127, 255}));
  connect(roomCycle.port_b, BT.fluPorVol[5]) annotation(
      Line(points = {{62, 82}, {40, 82}, {40, 70}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-180, 120}, {160, -60}})),
      Icon(graphics = {Rectangle(fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {1, 2}, extent = {{-63, 38}, {63, -38}}, textString = "%name", textStyle = {TextStyle.Bold})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end plant2;
  
  model controller2
    extends Modelica.Icons.Example;
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 273.15 + 20) annotation(
      Placement(visible = true, transformation(origin = {-84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback annotation(
      Placement(visible = true, transformation(origin = {-56, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.PI pi(T = 300, initType = Modelica.Blocks.Types.Init.InitialState, k = 5) annotation(
      Placement(visible = true, transformation(origin = {-28, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	ATW.plant2 plant2(redeclare package Medium = Medium, redeclare package MediumAir = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {2, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(const.y, feedback.u1) annotation(
      Line(points = {{-72, 18}, {-64, 18}}, color = {0, 0, 127}));
    connect(feedback.y, pi.u) annotation(
      Line(points = {{-46, 18}, {-40, 18}}, color = {0, 0, 127}));
  connect(pi.y, plant2.f) annotation(
      Line(points = {{-16, 18}, {-8, 18}}, color = {0, 0, 127}));
  connect(plant2.T, feedback.u2) annotation(
      Line(points = {{12, 18}, {20, 18}, {20, -6}, {-56, -6}, {-56, 10}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
      Diagram(coordinateSystem(extent = {{-100, 80}, {20, -20}})));
  end controller2;
  annotation(
    uses(Modelica(version = "4.0.0"), Buildings(version = "9.1.0")));
end ATW;