package A
  model FlowResistance_ex01
    replaceable package Medium = Media.Water.StandardWaterOnePhase;
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 10, offset = 10, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-90, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, T = 273.15 + 15, m_flow = 1, use_m_flow_in = true, nPorts = 1) annotation(
      Placement(visible = true, transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.05, length = 5, roughness = 4e-05) annotation(
      Placement(visible = true, transformation(origin = {10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 15 + 273.15, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {78, 10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = 273.15 + 15, V = 1e-3, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, p_start = 101325, use_T_start = true, use_portsData = false, nPorts = 2) annotation(
      Placement(visible = true, transformation(origin = {-20, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Medium, T_start = 273.15 + 15, V = 1e-3, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, p_start = 101325, use_HeatTransfer = false, use_T_start = true, use_portsData = false, nPorts = 2) annotation(
      Placement(visible = true, transformation(origin = {44, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.m_flow_in, ramp.y) annotation(
      Line(points = {{-60, 18}, {-72, 18}, {-72, 10}, {-78, 10}}, color = {0, 0, 127}));
    connect(boundary.ports[1], volume.ports[1]) annotation(
      Line(points = {{-40, 10}, {-20, 10}, {-20, 14}}, color = {0, 127, 255}));
    connect(volume.ports[2], pipe.port_a) annotation(
      Line(points = {{-20, 14}, {-18, 14}, {-18, 10}, {0, 10}}, color = {0, 127, 255}));
    connect(pipe.port_b, volume1.ports[1]) annotation(
      Line(points = {{20, 10}, {44, 10}, {44, 14}}, color = {0, 127, 255}));
    connect(volume1.ports[2], boundary1.ports[1]) annotation(
      Line(points = {{44, 14}, {46, 14}, {46, 10}, {68, 10}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
  end FlowResistance_ex01;

  import Modelica.Media;

  model FlowWithHeating_ex01
    //  extends Modelica.Icons.Example;
    //----------
    //replaceable package liquid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = liquid1
    //----------
    inner Modelica.Fluid.System system(T_ambient(displayUnit = "K") = 15 + 273.15, p_ambient(displayUnit = "Pa") = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 15 + 273.15, m_flow = 1, nPorts = 1, use_m_flow_in = false) annotation(
      Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 15 + 273.15, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {80, 10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp_heat(duration = 10, height = 100*1000, offset = 0, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1 annotation(
      Placement(visible = true, transformation(origin = {-40, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, V = Modelica.Constants.pi/4*pipe.diameter^2.0*pipe.length, nPorts = 2, use_HeatTransfer = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-10, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, diameter = 0.01, length = 0.5) annotation(
      Placement(visible = true, transformation(origin = {-40, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, diameter = 0.01, flowModel(show_Res = true), length = 0.5) annotation(
      Placement(visible = true, transformation(origin = {20, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, V = Modelica.Constants.pi/4*pipe1.diameter^2.0*pipe1.length, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {50, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ramp_heat.y, prescribedHeatFlow1.Q_flow) annotation(
      Line(points = {{-59, 50}, {-50, 50}}, color = {0, 0, 127}));
    connect(prescribedHeatFlow1.port, volume.heatPort) annotation(
      Line(points = {{-30, 50}, {-20, 50}, {-20, 20}}, color = {191, 0, 0}));
    connect(boundary.ports[1], pipe.port_a) annotation(
      Line(points = {{-60, 10}, {-50, 10}, {-50, 10}, {-50, 10}}, color = {0, 127, 255}));
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{-30, 10}, {-12, 10}, {-12, 10}, {-10, 10}}, color = {0, 127, 255}));
    connect(volume.ports[2], pipe1.port_a) annotation(
      Line(points = {{-10, 10}, {10, 10}, {10, 10}, {10, 10}}, color = {0, 127, 255}));
    connect(pipe1.port_b, volume1.ports[1]) annotation(
      Line(points = {{30, 10}, {46, 10}, {46, 10}, {50, 10}}, color = {0, 127, 255}));
    connect(volume1.ports[2], boundary1.ports[1]) annotation(
      Line(points = {{50, 10}, {70, 10}}, color = {0, 127, 255}, thickness = 0.5));
    annotation(
      experiment(StartTime = 0, StopTime = 50, Tolerance = 1e-06, Interval = 0.1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-100, -100}, {120, 100}})),
      __OpenModelica_commandLineOptions = "");
  end FlowWithHeating_ex01;

  model pumpingSystem_ex01
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid1
    //----------
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 288.15, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, diameter = 0.05, length = 0.5, nParallel = 1) annotation(
      Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {80, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 1000, offset = 1000, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start, use_N_in = true) annotation(
      Placement(visible = true, transformation(origin = {-40, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.RelativePressure relativePressure1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {-40, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {-10, 30}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Fluid.Sensors.MassFlowRate massFlowRate1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {20, 30}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Blocks.Math.Product d_flowPwr annotation(
      Placement(visible = true, transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Math.Gain gain_effPump(k = 1/0.9) annotation(
      Placement(visible = true, transformation(origin = {-30, -80}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(pump.port_b, relativePressure1.port_a) annotation(
      Line(points = {{-30, 30}, {-26, 30}, {-26, -20}, {-30, -20}}, color = {0, 127, 255}));
    connect(pump.port_a, relativePressure1.port_b) annotation(
      Line(points = {{-50, 30}, {-58, 30}, {-58, -20}, {-50, -20}}, color = {0, 127, 255}));
    connect(relativePressure1.p_rel, d_flowPwr.u2) annotation(
      Line(points = {{-40, -29}, {-40, -32}, {-36, -32}, {-36, -38}}, color = {0, 0, 127}));
    connect(boundary.ports[1], pump.port_a) annotation(
      Line(points = {{-70, 30}, {-50, 30}, {-50, 30}, {-50, 30}}, color = {0, 127, 255}));
    connect(ramp1.y, pump.N_in) annotation(
      Line(points = {{-59, 60}, {-40, 60}, {-40, 40}}, color = {0, 0, 127}));
    connect(pump.port_b, volumeFlowRate1.port_a) annotation(
      Line(points = {{-30, 30}, {-20, 30}, {-20, 30}, {-20, 30}}, color = {0, 127, 255}));
    connect(volumeFlowRate1.port_b, massFlowRate1.port_a) annotation(
      Line(points = {{0, 30}, {10, 30}, {10, 30}, {10, 30}}, color = {0, 127, 255}));
    connect(volumeFlowRate1.V_flow, d_flowPwr.u1) annotation(
      Line(points = {{-10, 19}, {-10, -30}, {-24, -30}, {-24, -38}}, color = {0, 0, 127}));
    connect(massFlowRate1.port_b, pipe.port_a) annotation(
      Line(points = {{30, 30}, {40, 30}, {40, 30}, {40, 30}}, color = {0, 127, 255}));
    connect(pipe.port_b, boundary1.ports[1]) annotation(
      Line(points = {{60, 30}, {70, 30}}, color = {0, 127, 255}));
    connect(gain_effPump.u, d_flowPwr.y) annotation(
      Line(points = {{-30, -68}, {-30, -68}, {-30, -60}, {-30, -60}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 40, Tolerance = 1e-06, Interval = 0.08),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      Diagram,
      __OpenModelica_commandLineOptions = "");
  end pumpingSystem_ex01;

  model pumpingSystem_ex02
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid1
    //----------
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 288.15, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-110, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, diameter = 0.05, length = 0.5, nParallel = 1) annotation(
      Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {90, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 1000, offset = 1000, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.UnitConversions.From_rpm from_rpm1 annotation(
      Placement(visible = true, transformation(origin = {-80, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.Pump pump1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start) annotation(
      Placement(visible = true, transformation(origin = {-30, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed1(phi(fixed = false)) annotation(
      Placement(visible = true, transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation(
      Placement(visible = true, transformation(origin = {-30, 20}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Fluid.Sensors.RelativePressure relativePressure1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {-30, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain gain_pumpEff(k = 0.9) annotation(
      Placement(visible = true, transformation(origin = {0, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Product d_flowPwrOutlet annotation(
      Placement(visible = true, transformation(origin = {-24, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate2(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {-70, -10}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Fluid.Sensors.MassFlowRate massFlowRate1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Blocks.Math.Product d_flowPwrInlet annotation(
      Placement(visible = true, transformation(origin = {-64, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(d_flowPwrOutlet.u1, volumeFlowRate1.V_flow) annotation(
      Line(points = {{-18, -88}, {-18, -78}, {0, -78}, {0, -22}}, color = {0, 0, 127}));
    connect(relativePressure1.port_a, pump1.port_b) annotation(
      Line(points = {{-20, -70}, {-14, -70}, {-14, -10}, {-20, -10}}, color = {0, 127, 255}));
    connect(pump1.port_a, relativePressure1.port_b) annotation(
      Line(points = {{-40, -10}, {-48, -10}, {-48, -70}, {-40, -70}}, color = {0, 127, 255}));
    connect(relativePressure1.p_rel, d_flowPwrOutlet.u2) annotation(
      Line(points = {{-30, -79}, {-30, -87}}, color = {0, 0, 127}));
    connect(relativePressure1.p_rel, d_flowPwrInlet.u1) annotation(
      Line(points = {{-30, -79}, {-30, -79}, {-30, -83}, {-58, -83}, {-58, -89}, {-58, -89}}, color = {0, 0, 127}));
    connect(volumeFlowRate2.V_flow, d_flowPwrInlet.u2) annotation(
      Line(points = {{-70, -20}, {-70, -88}}, color = {0, 0, 127}));
    connect(boundary.ports[1], volumeFlowRate2.port_a) annotation(
      Line(points = {{-100, -10}, {-80, -10}}, color = {0, 127, 255}));
    connect(pipe.port_b, boundary1.ports[1]) annotation(
      Line(points = {{70, -10}, {80, -10}}, color = {0, 127, 255}));
    connect(massFlowRate1.port_b, pipe.port_a) annotation(
      Line(points = {{40, -10}, {50, -10}}, color = {0, 127, 255}));
    connect(volumeFlowRate1.port_b, massFlowRate1.port_a) annotation(
      Line(points = {{10, -10}, {20, -10}, {20, -10}, {20, -10}}, color = {0, 127, 255}));
    connect(volumeFlowRate2.port_b, pump1.port_a) annotation(
      Line(points = {{-60, -10}, {-40, -10}, {-40, -10}, {-40, -10}}, color = {0, 127, 255}));
    connect(powerSensor1.power, gain_pumpEff.u) annotation(
      Line(points = {{-18, 28}, {-12, 28}, {-12, 28}, {-12, 28}}, color = {0, 0, 127}));
    connect(pump1.port_b, volumeFlowRate1.port_a) annotation(
      Line(points = {{-20, -10}, {-10, -10}}, color = {0, 127, 255}));
    connect(powerSensor1.flange_b, pump1.shaft) annotation(
      Line(points = {{-30, 10}, {-30, 0}}));
    connect(speed1.flange, powerSensor1.flange_a) annotation(
      Line(points = {{-40, 50}, {-30, 50}, {-30, 30}}));
    connect(from_rpm1.y, speed1.w_ref) annotation(
      Line(points = {{-69, 50}, {-63, 50}, {-63, 50}, {-63, 50}}, color = {0, 0, 127}));
    connect(ramp1.y, from_rpm1.u) annotation(
      Line(points = {{-99, 50}, {-93, 50}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 40, Tolerance = 1e-06, Interval = 0.08),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-140, -120}, {120, 100}})),
      __OpenModelica_commandLineOptions = "");
  end pumpingSystem_ex02;

  model pump_heater
    replaceable package fluid = Media.Water.StandardWaterOnePhase;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump1(redeclare package Medium = fluid, N_nominal = 1750, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.034, 0.04}, head_nominal = {39.0, 27.0, 22.8}), m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_b_start = 10*system.p_start, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower(V_flow_nominal = {0, 0.034, 0.04}, W_nominal = {5000, 14700, 17000}), use_N_in = true, use_powerCharacteristic = true) annotation(
      Placement(visible = true, transformation(origin = {0, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1000, offset = 1000, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-34, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = fluid, diameter = 0.05, length = 0.5, nParallel = 1) annotation(
      Placement(visible = true, transformation(origin = {32, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = fluid, T_start = 15 + 273.15, V = 0.1, nPorts = 2, p_start = 101325, use_HeatTransfer = false, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-48, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_ph boundary_start(redeclare package Medium = fluid, nPorts = 1, use_h_in = true, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-60, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = fluid, nPorts = 1, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {4, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101325) annotation(
      Placement(visible = true, transformation(origin = {-84, -64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {36, -84}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  equation
    connect(ramp.y, pump1.N_in) annotation(
      Line(points = {{-22, 50}, {0, 50}, {0, 16}}, color = {0, 0, 127}));
    connect(pump1.port_b, pipe.port_a) annotation(
      Line(points = {{10, 6}, {22, 6}}, color = {0, 127, 255}));
    connect(boundary_start.ports[1], volume.ports[1]) annotation(
      Line(points = {{-60, -16}, {-60, 6}, {-48, 6}}, color = {0, 127, 255}));
    connect(volume.ports[2], pump1.port_a) annotation(
      Line(points = {{-48, 6}, {-10, 6}}, color = {0, 127, 255}));
    connect(const_p_tank.y, boundary_start.p_in) annotation(
      Line(points = {{-73, -64}, {-68, -64}, {-68, -38}}, color = {0, 0, 127}));
    connect(boundary.p_in, const_p_tank.y) annotation(
      Line(points = {{-8, -64}, {-73, -64}}, color = {0, 0, 127}));
    connect(pipe.port_b, specificEnthalpy.port) annotation(
      Line(points = {{42, 6}, {64, 6}, {64, -74}, {36, -74}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port, boundary.ports[1]) annotation(
      Line(points = {{36, -74}, {25, -74}, {25, -72}, {14, -72}}, color = {0, 127, 255}));
    connect(specificEnthalpy.h_out, boundary_start.h_in) annotation(
      Line(points = {{25, -84}, {-64, -84}, {-64, -38}}, color = {0, 0, 127}));
    annotation(
      Diagram(graphics = {Text(origin = {20, 60}, extent = {{-6, 4}, {6, -4}}, textString = "text")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end pump_heater;

  model pumpheater2
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid1
    //----------
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 288.15, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, diameter = 0.05, length = 0.5, nParallel = 1) annotation(
      Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, nPorts = 1, p = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {80, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 1000, offset = 1000, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start, use_N_in = true) annotation(
      Placement(visible = true, transformation(origin = {-40, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.ports[1], pump.port_a) annotation(
      Line(points = {{-70, 30}, {-50, 30}, {-50, 30}, {-50, 30}}, color = {0, 127, 255}));
    connect(ramp1.y, pump.N_in) annotation(
      Line(points = {{-59, 60}, {-40, 60}, {-40, 40}}, color = {0, 0, 127}));
    connect(pipe.port_b, boundary1.ports[1]) annotation(
      Line(points = {{60, 30}, {70, 30}}, color = {0, 127, 255}));
    connect(pump.port_b, pipe.port_a) annotation(
      Line(points = {{-30, 30}, {40, 30}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 40, Tolerance = 1e-06, Interval = 0.08),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
      Diagram,
      __OpenModelica_commandLineOptions = "");
  end pumpheater2;

  model CoolingSystem_ex01
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Glycol47;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Essotherm650;
    //redeclare package Medium = fluid1
    //-----
    //replaceable package fluid2 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid2
    //----------
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp_pump_N(duration = 10, height = 0, offset = 500, startTime = 30) annotation(
      Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.UnitConversions.From_rpm from_rpm1 annotation(
      Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.Pump pump1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed1(phi(fixed = false)) annotation(
      Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation(
      Placement(visible = true, transformation(origin = {-50, 18}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe cooler_hside(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 10, use_HeatTransfer = true) annotation(
      Placement(visible = true, transformation(origin = {120, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe pipe(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 2, use_HeatTransfer = false) annotation(
      Placement(visible = true, transformation(origin = {184, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary2(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T = 15 + 273.15, m_flow = 5, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {10, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe cooler_cside(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T_start = 15 + 273.15, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = cooler_hside.nNodes, use_HeatTransfer = true) annotation(
      Placement(visible = true, transformation(origin = {120, 76}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary3(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, nPorts = 1) annotation(
      Placement(visible = true, transformation(origin = {220, 110}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection1[cooler_hside.nNodes] annotation(
      Placement(visible = true, transformation(origin = {120, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1[cooler_hside.nNodes](each C = 10) annotation(
      Placement(visible = true, transformation(origin = {130, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection2[cooler_hside.nNodes] annotation(
      Placement(visible = true, transformation(origin = {120, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Routing.Replicator replicator1(nout = cooler_hside.nNodes) annotation(
      Placement(visible = true, transformation(origin = {90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Routing.Replicator replicator2(nout = cooler_cside.nNodes) annotation(
      Placement(visible = true, transformation(origin = {90, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 5000) annotation(
      Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, V = 1*0.001, nPorts = 2, use_HeatTransfer = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {210, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1 annotation(
      Placement(visible = true, transformation(origin = {210, -38}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Ramp ramp_heat_generation(duration = 10, height = 10*1000, offset = 100*1000, startTime = 100) annotation(
      Placement(visible = true, transformation(origin = {198, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe pipe1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 2, use_HeatTransfer = false) annotation(
      Placement(visible = true, transformation(origin = {170, -110}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_ph boundary(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_h_in = true, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, T_start = 15 + 273.15, V = 20*0.001, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {70, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {110, -120}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase, T_start = 15 + 273.15, V = 1*0.001, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {180, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp_m_flow_coolant(duration = 10, height = 0.2, offset = 2, startTime = 160) annotation(
      Placement(visible = true, transformation(origin = {-30, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {4, -60}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-10, 27}, extent = {{-12, -9}, {12, 9}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue5(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {90, 100}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue6(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {160, 100}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue7(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {156, -30}, extent = {{12, -8}, {-12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue8(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {219, -109}, extent = {{-13, -7}, {13, 7}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue9(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {144, -92}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-140, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {-80, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {-20, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature2(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {158, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature3(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {200, -100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature4(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {130, -110}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature5(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {60, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature6(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {130, 100}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor[cooler_hside.nNodes] annotation(
      Placement(visible = true, transformation(origin = {120, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Math.Sum sum1(nin = cooler_hside.nNodes) annotation(
      Placement(visible = true, transformation(origin = {150, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue4(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {2, -76}, extent = {{12, -7}, {-12, 7}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue3(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {32, -28}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue2(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-22, -32}, extent = {{12, -8}, {-12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue1(significantDigits = 5, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-64, -28}, extent = {{-26, -8}, {26, 8}}, rotation = 0)));
  equation
    connect(volume1.ports[1], boundary3.ports[1]) annotation(
      Line(points = {{180, 110}, {210, 110}}, color = {0, 127, 255}, thickness = 0.5));
    connect(const.y, replicator1.u) annotation(
      Line(points = {{71, -10}, {73, -10}, {73, -30}, {78, -30}}, color = {0, 0, 127}));
    connect(const.y, replicator2.u) annotation(
      Line(points = {{71, -10}, {73, -10}, {73, 18}, {78, 18}}, color = {0, 0, 127}));
    connect(replicator2.y, convection2.Gc) annotation(
      Line(points = {{101, 18}, {110, 18}}, color = {0, 0, 127}, thickness = 0.5));
    connect(heatCapacitor1.port, convection2.solid) annotation(
      Line(points = {{120, -6}, {120, 8}}, color = {191, 0, 0}, thickness = 0.5));
    connect(replicator1.y, convection1.Gc) annotation(
      Line(points = {{101, -30}, {110, -30}}, color = {0, 0, 127}, thickness = 0.5));
    connect(heatCapacitor1.port, convection1.solid) annotation(
      Line(points = {{120, -6}, {120, -20}}, color = {191, 0, 0}, thickness = 0.5));
    connect(convection1.fluid, cooler_hside.heatPorts) annotation(
      Line(points = {{120, -40}, {120, -46}}, color = {191, 0, 0}, thickness = 0.5));
    connect(speed1.flange, powerSensor1.flange_a) annotation(
      Line(points = {{-60, 50}, {-50, 50}, {-50, 28}}));
    connect(powerSensor1.flange_b, pump1.shaft) annotation(
      Line(points = {{-50, 8}, {-50, -40}}));
    connect(from_rpm1.y, speed1.w_ref) annotation(
      Line(points = {{-89, 50}, {-83, 50}, {-83, 50}, {-83, 50}}, color = {0, 0, 127}));
    connect(ramp_pump_N.y, from_rpm1.u) annotation(
      Line(points = {{-119, 70}, {-116, 70}, {-116, 50}, {-113, 50}}, color = {0, 0, 127}));
    connect(ramp_heat_generation.y, prescribedHeatFlow1.Q_flow) annotation(
      Line(points = {{209, -10}, {210, -10}, {210, -28}}, color = {0, 0, 127}));
    connect(prescribedHeatFlow1.port, volume.heatPort) annotation(
      Line(points = {{210, -48}, {210, -60}}, color = {191, 0, 0}));
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{194, -50}, {200, -50}, {200, -70}}, color = {0, 127, 255}));
    connect(boundary.ports[1], tank.ports[1]) annotation(
      Line(points = {{-110, -70}, {-110, -50}}, color = {0, 127, 255}, thickness = 0.5));
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{80, -100}, {95, -100}, {95, -110}, {110, -110}}, color = {0, 127, 255}));
    connect(volumeFlowRate1.port_b, cooler_hside.port_a) annotation(
      Line(points = {{40, -50}, {110, -50}}, color = {0, 127, 255}));
    connect(ramp_m_flow_coolant.y, boundary2.m_flow_in) annotation(
      Line(points = {{-19, 120}, {-17, 120}, {-17, 118}, {-1, 118}}, color = {0, 0, 127}));
    connect(pressure.port, volumeFlowRate1.port_a) annotation(
      Line(points = {{4, -50}, {20, -50}}, color = {0, 127, 255}));
    connect(powerSensor1.power, realValue.numberPort) annotation(
      Line(points = {{-39, 26}, {-25.75, 26}, {-25.75, 27}, {-24, 27}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary1.p_in) annotation(
      Line(points = {{-129, -110}, {-33.5, -110}, {-33.5, -92}, {58, -92}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary.p_in) annotation(
      Line(points = {{-128, -110}, {-118, -110}, {-118, -92}}, color = {0, 0, 127}));
    connect(specificEnthalpy.h_out, boundary.h_in) annotation(
      Line(points = {{100, -120}, {-114, -120}, {-114, -92}}, color = {0, 0, 127}));
    connect(tank.ports[2], temperature.port_a) annotation(
      Line(points = {{-110, -50}, {-90, -50}}, color = {0, 127, 255}));
    connect(temperature.port_b, pump1.port_a) annotation(
      Line(points = {{-70, -50}, {-60, -50}}, color = {0, 127, 255}));
    connect(pump1.port_b, temperature1.port_a) annotation(
      Line(points = {{-40, -50}, {-30, -50}}));
    connect(temperature1.port_b, pressure.port) annotation(
      Line(points = {{-10, -50}, {4, -50}}, color = {0, 127, 255}));
    connect(cooler_hside.port_b, temperature2.port_a) annotation(
      Line(points = {{130, -50}, {148, -50}}, color = {0, 127, 255}));
    connect(temperature2.port_b, pipe.port_a) annotation(
      Line(points = {{168, -50}, {174, -50}}, color = {0, 127, 255}));
    connect(temperature2.T, realValue7.numberPort) annotation(
      Line(points = {{158, -38}, {170, -38}, {170, -30}}, color = {0, 0, 127}));
    connect(pipe1.port_a, temperature3.port_b) annotation(
      Line(points = {{180, -110}, {200, -110}}, color = {0, 127, 255}));
    connect(volume.ports[2], temperature3.port_a) annotation(
      Line(points = {{200, -70}, {200, -90}}, color = {0, 127, 255}));
    connect(temperature3.T, realValue8.numberPort) annotation(
      Line(points = {{212, -100}, {212, -106.5}, {204, -106.5}, {204, -109}}, color = {0, 0, 127}));
    connect(temperature4.port_a, pipe1.port_b) annotation(
      Line(points = {{140, -110}, {160, -110}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port, temperature4.port_b) annotation(
      Line(points = {{110, -110}, {120, -110}}, color = {0, 127, 255}));
    connect(temperature4.T, realValue9.numberPort) annotation(
      Line(points = {{130, -98}, {130, -92}}, color = {0, 0, 127}));
    connect(boundary2.ports[1], temperature5.port_a) annotation(
      Line(points = {{20, 110}, {60, 110}}, color = {0, 127, 255}));
    connect(temperature5.T, realValue5.numberPort) annotation(
      Line(points = {{71, 100}, {75, 100}}, color = {0, 0, 127}));
    connect(temperature5.port_b, cooler_cside.port_a) annotation(
      Line(points = {{60, 90}, {60, 76}, {110, 76}}, color = {0, 127, 255}));
    connect(cooler_cside.port_b, temperature6.port_a) annotation(
      Line(points = {{130, 76}, {130, 90}}, color = {0, 127, 255}));
    connect(temperature6.port_b, volume1.ports[2]) annotation(
      Line(points = {{130, 110}, {180, 110}}, color = {0, 127, 255}));
    connect(temperature6.T, realValue6.numberPort) annotation(
      Line(points = {{141, 100}, {145, 100}}, color = {0, 0, 127}));
    connect(heatFlowSensor.port_b, cooler_cside.heatPorts) annotation(
      Line(points = {{120, 60}, {120, 72}}, color = {191, 0, 0}, thickness = 0.5));
    connect(heatFlowSensor.port_a, convection2.fluid) annotation(
      Line(points = {{120, 40}, {120, 28}}, color = {191, 0, 0}, thickness = 0.5));
    connect(heatFlowSensor.Q_flow, sum1.u) annotation(
      Line(points = {{130, 50}, {138, 50}}, color = {0, 0, 127}, thickness = 0.5));
    connect(pressure.p, realValue4.numberPort) annotation(
      Line(points = {{15, -60}, {15, -76}}, color = {0, 0, 127}));
    connect(volumeFlowRate1.V_flow, realValue3.numberPort) annotation(
      Line(points = {{30, -39}, {30, -34}, {18, -34}, {18, -28}}, color = {0, 0, 127}));
    connect(temperature1.T, realValue2.numberPort) annotation(
      Line(points = {{-20, -38}, {-8, -38}, {-8, -32}}, color = {0, 0, 127}));
    connect(temperature.T, realValue1.numberPort) annotation(
      Line(points = {{-80, -38}, {-80, -33}, {-94, -33}, {-94, -28}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 260, Tolerance = 1e-06, Interval = 0.05),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, -140}, {240, 140}}, initialScale = 0.1), graphics = {Rectangle(origin = {110, 28}, extent = {{-62, 96}, {40, -101}}), Text(origin = {94, 128}, extent = {{-34, 4}, {34, -4}}, textString = "Heat Exchanger"), Text(origin = {2, 96}, extent = {{-34, 4}, {34, -4}}, textString = "coolant flow line")}),
      __OpenModelica_commandLineOptions = "");
  end CoolingSystem_ex01;

  model CoolingSystem_ex01_test
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Glycol47;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Essotherm650;
    //redeclare package Medium = fluid1
    //-----
    //replaceable package fluid2 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid2
    //----------
    Modelica.Blocks.Sources.Ramp ramp_pump_N(duration = 10, height = 0, offset = 500, startTime = 30) annotation(
      Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.UnitConversions.From_rpm from_rpm1 annotation(
      Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.Pump pump1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed1(phi(fixed = false)) annotation(
      Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation(
      Placement(visible = true, transformation(origin = {-50, 18}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Fluid.Sensors.VolumeFlowRate volumeFlowRate1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe pipe1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 2, use_HeatTransfer = false) annotation(
      Placement(visible = true, transformation(origin = {170, -110}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_ph boundary(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_h_in = true, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, T_start = 15 + 273.15, V = 20*0.001, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {70, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {110, -120}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase) annotation(
      Placement(visible = true, transformation(origin = {4, -60}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-10, 27}, extent = {{-12, -9}, {12, 9}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue9(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {144, -92}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-140, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {-80, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {-20, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature4(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {130, -110}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue4(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {2, -76}, extent = {{12, -7}, {-12, 7}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue3(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {32, -28}, extent = {{-12, -8}, {12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue2(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-22, -32}, extent = {{12, -8}, {-12, 8}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue1(significantDigits = 5, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-64, -28}, extent = {{-26, -8}, {26, 8}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {212, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(speed1.flange, powerSensor1.flange_a) annotation(
      Line(points = {{-60, 50}, {-50, 50}, {-50, 28}}));
    connect(powerSensor1.flange_b, pump1.shaft) annotation(
      Line(points = {{-50, 8}, {-50, -40}}));
    connect(from_rpm1.y, speed1.w_ref) annotation(
      Line(points = {{-89, 50}, {-83, 50}, {-83, 50}, {-83, 50}}, color = {0, 0, 127}));
    connect(ramp_pump_N.y, from_rpm1.u) annotation(
      Line(points = {{-119, 70}, {-116, 70}, {-116, 50}, {-113, 50}}, color = {0, 0, 127}));
    connect(boundary.ports[1], tank.ports[1]) annotation(
      Line(points = {{-110, -70}, {-110, -50}}, color = {0, 127, 255}, thickness = 0.5));
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{80, -100}, {95, -100}, {95, -110}, {110, -110}}, color = {0, 127, 255}));
    connect(pressure.port, volumeFlowRate1.port_a) annotation(
      Line(points = {{4, -50}, {20, -50}}, color = {0, 127, 255}));
    connect(powerSensor1.power, realValue.numberPort) annotation(
      Line(points = {{-39, 26}, {-25.75, 26}, {-25.75, 27}, {-24, 27}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary1.p_in) annotation(
      Line(points = {{-129, -110}, {-33.5, -110}, {-33.5, -92}, {58, -92}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary.p_in) annotation(
      Line(points = {{-128, -110}, {-118, -110}, {-118, -92}}, color = {0, 0, 127}));
    connect(specificEnthalpy.h_out, boundary.h_in) annotation(
      Line(points = {{100, -120}, {-114, -120}, {-114, -92}}, color = {0, 0, 127}));
    connect(tank.ports[2], temperature.port_a) annotation(
      Line(points = {{-110, -50}, {-90, -50}}, color = {0, 127, 255}));
    connect(temperature.port_b, pump1.port_a) annotation(
      Line(points = {{-70, -50}, {-60, -50}}, color = {0, 127, 255}));
    connect(pump1.port_b, temperature1.port_a) annotation(
      Line(points = {{-40, -50}, {-30, -50}}));
    connect(temperature1.port_b, pressure.port) annotation(
      Line(points = {{-10, -50}, {4, -50}}, color = {0, 127, 255}));
    connect(temperature4.port_a, pipe1.port_b) annotation(
      Line(points = {{140, -110}, {160, -110}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port, temperature4.port_b) annotation(
      Line(points = {{110, -110}, {120, -110}}, color = {0, 127, 255}));
    connect(temperature4.T, realValue9.numberPort) annotation(
      Line(points = {{130, -98}, {130, -92}}, color = {0, 0, 127}));
    connect(pressure.p, realValue4.numberPort) annotation(
      Line(points = {{15, -60}, {15, -76}}, color = {0, 0, 127}));
    connect(volumeFlowRate1.V_flow, realValue3.numberPort) annotation(
      Line(points = {{30, -39}, {30, -34}, {18, -34}, {18, -28}}, color = {0, 0, 127}));
    connect(temperature1.T, realValue2.numberPort) annotation(
      Line(points = {{-20, -38}, {-8, -38}, {-8, -32}}, color = {0, 0, 127}));
    connect(temperature.T, realValue1.numberPort) annotation(
      Line(points = {{-80, -38}, {-80, -33}, {-94, -33}, {-94, -28}}, color = {0, 0, 127}));
    connect(volumeFlowRate1.port_b, pipe1.port_a) annotation(
      Line(points = {{40, -50}, {188, -50}, {188, -110}, {180, -110}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 260, Tolerance = 1e-06, Interval = 0.05),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, -140}, {240, 140}}, initialScale = 0.1)),
      __OpenModelica_commandLineOptions = "");
  end CoolingSystem_ex01_test;

  model CoolingSystem_ex01_test1
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Glycol47;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Essotherm650;
    //redeclare package Medium = fluid1
    //-----
    //replaceable package fluid2 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid2
    //----------
    Modelica.Blocks.Sources.Ramp ramp_pump_N(duration = 10, height = 0, offset = 500, startTime = 30) annotation(
      Placement(visible = true, transformation(origin = {-130, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.UnitConversions.From_rpm from_rpm1 annotation(
      Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Machines.Pump pump1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, N_nominal = 1000, V(displayUnit = "l") = 0.001, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start) annotation(
      Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed1(phi(fixed = false)) annotation(
      Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation(
      Placement(visible = true, transformation(origin = {-50, 18}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Fluid.Sources.Boundary_ph boundary(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_h_in = true, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, T_start = 15 + 273.15, V = 20*0.001, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-14, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {24, -114}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Blocks.Interaction.Show.RealValue realValue(significantDigits = 4, use_numberPort = true) annotation(
      Placement(visible = true, transformation(origin = {-10, 27}, extent = {{-12, -9}, {12, 9}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-140, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {212, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2) annotation(
      Placement(visible = true, transformation(origin = {36, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe pipe1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 2, use_HeatTransfer = false) annotation(
      Placement(visible = true, transformation(origin = {-18, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
  equation
    connect(speed1.flange, powerSensor1.flange_a) annotation(
      Line(points = {{-60, 50}, {-50, 50}, {-50, 28}}));
    connect(powerSensor1.flange_b, pump1.shaft) annotation(
      Line(points = {{-50, 8}, {-50, -40}}));
    connect(from_rpm1.y, speed1.w_ref) annotation(
      Line(points = {{-89, 50}, {-83, 50}, {-83, 50}, {-83, 50}}, color = {0, 0, 127}));
    connect(ramp_pump_N.y, from_rpm1.u) annotation(
      Line(points = {{-119, 70}, {-116, 70}, {-116, 50}, {-113, 50}}, color = {0, 0, 127}));
    connect(boundary.ports[1], tank.ports[1]) annotation(
      Line(points = {{-110, -70}, {-110, -50}}, color = {0, 127, 255}, thickness = 0.5));
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{-4, -94}, {10, -94}, {10, -104}, {24, -104}}, color = {0, 127, 255}));
    connect(powerSensor1.power, realValue.numberPort) annotation(
      Line(points = {{-39, 26}, {-25.75, 26}, {-25.75, 27}, {-24, 27}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary1.p_in) annotation(
      Line(points = {{-129, -110}, {-33.5, -110}, {-33.5, -86}, {-26, -86}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary.p_in) annotation(
      Line(points = {{-128, -110}, {-118, -110}, {-118, -92}}, color = {0, 0, 127}));
    connect(tank.ports[2], pump1.port_a) annotation(
      Line(points = {{-110, -50}, {-60, -50}}, color = {0, 127, 255}));
    connect(boundary.h_in, specificEnthalpy.h_out) annotation(
      Line(points = {{-114, -92}, {-114, -114}, {14, -114}}, color = {0, 0, 127}));
    connect(pipe.port_b, specificEnthalpy.port) annotation(
      Line(points = {{46, -50}, {52, -50}, {52, -104}, {24, -104}}, color = {0, 127, 255}));
    connect(pipe1.port_b, pipe.port_a) annotation(
      Line(points = {{-8, -50}, {26, -50}}, color = {0, 127, 255}));
    connect(pipe1.port_a, pump1.port_b) annotation(
      Line(points = {{-28, -50}, {-40, -50}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 260, Tolerance = 1e-06, Interval = 0.05),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, -140}, {240, 140}}, initialScale = 0.1)),
      __OpenModelica_commandLineOptions = "");
  end CoolingSystem_ex01_test1;

  model CoolingSystem_ex01_test2
    extends Modelica.Icons.Example;
    //----------
    //replaceable package fluid1 = Modelica.Media.Water.StandardWaterOnePhase;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Glycol47;
    //replaceable package fluid1 = Modelica.Media.Incompressible.Examples.Essotherm650;
    //redeclare package Medium = fluid1
    //-----
    //replaceable package fluid2 = Modelica.Media.Water.StandardWaterOnePhase;
    //redeclare package Medium = fluid2
    //----------
    Modelica.Fluid.Sources.Boundary_ph boundary(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_h_in = true, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, T_start = 15 + 273.15, V = 1, nPorts = 2, use_HeatTransfer = false, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, nPorts = 1, p = 101.325*1000, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-20, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47) annotation(
      Placement(visible = true, transformation(origin = {24, -114}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-140, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {36, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.DynamicPipe pipe1(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, diameter = 0.02, length = 2, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 2, nParallel = 1, use_HeatTransfer = false) annotation(
      Placement(visible = true, transformation(origin = {-18, -50}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Modelica.Media.Incompressible.Examples.Glycol47, N_nominal = 1000, V = 1, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 10, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, p_b_start = 10*system.p_start, use_N_in = true) annotation(
      Placement(visible = true, transformation(origin = {-54, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1000, offset = 1000, startTime = 10) annotation(
      Placement(visible = true, transformation(origin = {-80, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.ports[1], tank.ports[1]) annotation(
      Line(points = {{-110, -70}, {-110, -50}}, color = {0, 127, 255}, thickness = 0.5));
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{-10, -92}, {-21, -92}, {-21, -104}, {24, -104}}, color = {0, 127, 255}));
    connect(const_p_tank.y, boundary1.p_in) annotation(
      Line(points = {{-129, -110}, {-33.5, -110}, {-33.5, -84}, {-32, -84}}, color = {0, 0, 127}));
    connect(const_p_tank.y, boundary.p_in) annotation(
      Line(points = {{-128, -110}, {-118, -110}, {-118, -92}}, color = {0, 0, 127}));
    connect(boundary.h_in, specificEnthalpy.h_out) annotation(
      Line(points = {{-114, -92}, {-114, -114}, {14, -114}}, color = {0, 0, 127}));
    connect(ramp.y, pump.N_in) annotation(
      Line(points = {{-69, -20}, {-55, -20}, {-55, -40}}, color = {0, 0, 127}));
    connect(tank.ports[2], pump.port_a) annotation(
      Line(points = {{-110, -50}, {-64, -50}}, color = {0, 127, 255}));
    connect(pump.port_b, pipe1.port_a) annotation(
      Line(points = {{-44, -50}, {-28, -50}}, color = {0, 127, 255}));
    connect(pipe1.port_b, specificEnthalpy.port) annotation(
      Line(points = {{-8, -50}, {36, -50}, {36, -104}, {24, -104}}));
    annotation(
      experiment(StartTime = 0, StopTime = 260, Tolerance = 1e-06, Interval = 0.05),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, 0}, {60, -140}}, initialScale = 0.1)),
      __OpenModelica_commandLineOptions = "");
  end CoolingSystem_ex01_test2;

  model CoolingSystem_ex01_test3
    replaceable package fluid = Media.Water.StandardWaterOnePhase;
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = fluid, T_start = 50 + 273.15, V = 10, nPorts = 2, use_HeatTransfer = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-36, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = fluid, nPorts = 1, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {-14, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {24, -114}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Blocks.Sources.Constant const_p_tank(k = 101.325*1000) annotation(
      Placement(visible = true, transformation(origin = {-80, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {212, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = fluid, diameter = 0.02, length = 2, roughness = 0.005) annotation(
      Placement(visible = true, transformation(origin = {34, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_h boundary(redeclare package Medium = fluid, nPorts = 1, use_h_in = true, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-80, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = 1, offset = 0, startTime = 1) annotation(
      Placement(visible = true, transformation(origin = {-138, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
      Placement(visible = true, transformation(origin = {-82, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {-18, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue annotation(
      Placement(visible = true, transformation(origin = {28, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature1(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {74, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature2(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {-48, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue1 annotation(
      Placement(visible = true, transformation(origin = {-8, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue2 annotation(
      Placement(visible = true, transformation(origin = {120, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 10) annotation(
      Placement(visible = true, transformation(origin = {-138, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{-4, -94}, {10, -94}, {10, -104}, {24, -104}}, color = {0, 127, 255}));
    connect(const_p_tank.y, boundary1.p_in) annotation(
      Line(points = {{-69, -86}, {-26, -86}}, color = {0, 0, 127}));
    connect(specificEnthalpy.h_out, boundary.h_in) annotation(
      Line(points = {{14, -114}, {-108, -114}, {-108, -44}, {-92, -44}}, color = {0, 0, 127}));
    connect(ramp.y, boundary.m_flow_in) annotation(
      Line(points = {{-126, -22}, {-102, -22}, {-102, -40}, {-90, -40}}, color = {0, 0, 127}));
    connect(prescribedHeatFlow.port, tank.heatPort) annotation(
      Line(points = {{-72, -8}, {-62, -8}, {-62, -38}, {-46, -38}}, color = {191, 0, 0}));
    connect(temperature.T, realValue.numberPort) annotation(
      Line(points = {{-11, -30}, {16.5, -30}, {16.5, -24}}, color = {0, 0, 127}));
    connect(tank.ports[1], temperature.port) annotation(
      Line(points = {{-36, -48}, {-18, -48}, {-18, -40}}, color = {0, 127, 255}));
    connect(temperature.port, pipe.port_a) annotation(
      Line(points = {{-18, -40}, {-18, -48}, {24, -48}}, color = {0, 127, 255}));
    connect(pipe.port_b, temperature1.port) annotation(
      Line(points = {{44, -48}, {74, -48}, {74, -78}}, color = {0, 127, 255}));
    connect(temperature1.port, specificEnthalpy.port) annotation(
      Line(points = {{74, -78}, {82, -78}, {82, -104}, {24, -104}}, color = {0, 127, 255}));
    connect(boundary.ports[1], temperature2.port) annotation(
      Line(points = {{-70, -48}, {-64, -48}, {-64, -80}, {-48, -80}}, color = {0, 127, 255}));
    connect(temperature2.port, tank.ports[2]) annotation(
      Line(points = {{-48, -80}, {-30, -80}, {-30, -52}, {-50, -52}, {-50, -48}, {-36, -48}}, color = {0, 127, 255}));
    connect(temperature2.T, realValue1.numberPort) annotation(
      Line(points = {{-40, -70}, {-22, -70}, {-22, -68}, {-20, -68}}, color = {0, 0, 127}));
    connect(temperature1.T, realValue2.numberPort) annotation(
      Line(points = {{82, -68}, {98.25, -68}, {98.25, -76}, {108.5, -76}}, color = {0, 0, 127}));
    connect(const.y, prescribedHeatFlow.Q_flow) annotation(
      Line(points = {{-126, 12}, {-98, 12}, {-98, -8}, {-92, -8}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 50, Tolerance = 1e-06, Interval = 0.1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, -140}, {240, 140}}, initialScale = 0.1)),
      __OpenModelica_commandLineOptions = "");
  end CoolingSystem_ex01_test3;

  model heatExchanger
    replaceable package fluid = Media.Water.StandardWaterOnePhase;
    Modelica.Fluid.Vessels.ClosedVolume tank(redeclare package Medium = fluid, T_start = 273.15 + 20, V = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nPorts = 2, use_HeatTransfer = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {20, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = fluid, T = 273.15 + 20, nPorts = 1, use_p_in = false) annotation(
      Placement(visible = true, transformation(origin = {-14, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {212, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_h boundary(redeclare package Medium = fluid, nPorts = 1, use_h_in = false, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-80, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
      Placement(visible = true, transformation(origin = {-6, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {50, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue annotation(
      Placement(visible = true, transformation(origin = {80, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.Temperature temperature2(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {-42, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Interaction.Show.RealValue realValue1 annotation(
      Placement(visible = true, transformation(origin = {-16, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 1000) annotation(
      Placement(visible = true, transformation(origin = {-138, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const1(k = 1) annotation(
      Placement(visible = true, transformation(origin = {-136, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sensors.SpecificEnthalpy specificEnthalpy(redeclare package Medium = fluid) annotation(
      Placement(visible = true, transformation(origin = {26, -104}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  equation
    connect(prescribedHeatFlow.port, tank.heatPort) annotation(
      Line(points = {{4, -14}, {4, -38}, {10, -38}}, color = {191, 0, 0}));
    connect(temperature.T, realValue.numberPort) annotation(
      Line(points = {{57, -38}, {66.75, -38}, {66.75, -42}, {68.5, -42}}, color = {0, 0, 127}));
    connect(tank.ports[1], temperature.port) annotation(
      Line(points = {{20, -48}, {50, -48}}, color = {0, 127, 255}));
    connect(boundary.ports[1], temperature2.port) annotation(
      Line(points = {{-70, -48}, {-42, -48}}, color = {0, 127, 255}));
    connect(temperature2.port, tank.ports[2]) annotation(
      Line(points = {{-42, -48}, {20, -48}}, color = {0, 127, 255}));
    connect(temperature2.T, realValue1.numberPort) annotation(
      Line(points = {{-35, -38}, {-27.5, -38}}, color = {0, 0, 127}));
    connect(const.y, prescribedHeatFlow.Q_flow) annotation(
      Line(points = {{-126, 12}, {-98, 12}, {-98, -14}, {-16, -14}}, color = {0, 0, 127}));
    connect(const1.y, boundary.m_flow_in) annotation(
      Line(points = {{-125, -34}, {-107.5, -34}, {-107.5, -40}, {-90, -40}}, color = {0, 0, 127}));
    connect(boundary1.ports[1], specificEnthalpy.port) annotation(
      Line(points = {{-4, -94}, {26, -94}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port, temperature.port) annotation(
      Line(points = {{26, -94}, {60, -94}, {60, -48}, {50, -48}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 50, Tolerance = 1e-06, Interval = 0.1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
      Diagram(coordinateSystem(extent = {{-160, -140}, {240, 140}}, initialScale = 0.1)),
      __OpenModelica_commandLineOptions = "");
  end heatExchanger;

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
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor[pipe.nNodes](each C = pipe_rho*pi*(pipe_diameterOuter^2 - pipe_diameterInner^2)/4*pipe_length/pipe.nNodes*pipe_c, each T(start = 273.15 + 50)) annotation(
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
    A.pipe_ pipe_(redeclare package Medium = Medium) annotation(
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
    A.HEX_HP hex_hp(K = 2, redeclare package Medium = Medium, T = 100) annotation(
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
  	replaceable package Medium = Media.Air.DryAirNasa;
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
  
  	// def
  	Modelica.Fluid.Vessels.ClosedVolume air(redeclare package Medium = Medium, T_start = room_initialTemp, V = room_volume, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 0, use_HeatTransfer = true, use_T_start = true, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation(
      Placement(visible = true, transformation(origin = {-38, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Components.Convection convection_room annotation(
      Placement(visible = true, transformation(origin = {-2, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = wall_volume*wall_density*wall_specHeatCap, T(fixed = true, start = wall_initialTemp)) annotation(
      Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Components.Convection convection_amb annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T = T_amb)  annotation(
      Placement(visible = true, transformation(origin = {130, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  	Modelica.Blocks.Sources.Constant const_room(k = air_h*room_area) annotation(
      Placement(visible = true, transformation(origin = {-26, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = wall_TC*room_area*2/wall_thick)  annotation(
  		Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = wall_TC*room_area*2/wall_thick)  annotation(
  		Placement(visible = true, transformation(origin = {74, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  
  	// Interfaces
  	Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_rad annotation(
      Placement(visible = true, transformation(origin = {-100, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_conv annotation(
      Placement(visible = true, transformation(origin = {-100, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Blocks.Interfaces.RealOutput T annotation(
      Placement(visible = true, transformation(origin = {146, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
  	connect(port_rad, air.heatPort) annotation(
      Line(points = {{-100, 28}, {-62, 28}, {-62, 0}, {-52, 0}}, color = {191, 0, 0}));
  	connect(port_conv, air.heatPort) annotation(
      Line(points = {{-100, -18}, {-62, -18}, {-62, 0}, {-52, 0}}, color = {191, 0, 0}));
  	connect(air.heatPort, temperatureSensor.port) annotation(
      Line(points = {{-52, 0}, {-56, 0}, {-56, -42}, {-48, -42}}, color = {191, 0, 0}));
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
  		Line(points = {{-26, -42}, {146, -42}}, color = {0, 0, 127}));
  	annotation(
  		Icon(graphics = {Rectangle(fillColor = {170, 170, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -3}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Line(origin = {19.9746, 97.1102}, points = {{-62.0003, -38.0003}, {5.99968, 25.9997}, {-4.00032, -28.0003}, {61.9997, 37.9997}}, color = {255, 0, 0}, thickness = 0.75, arrow = {Arrow.None, Arrow.Open}, arrowSize = 15)}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
  		Diagram(coordinateSystem(extent = {{-120, 80}, {160, -60}})));
  end room;
  
  model test_room
  	replaceable package Medium = Media.Air.DryAirNasa;
  	inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1000, offset = 0, startTime = 5)  annotation(
      Placement(visible = true, transformation(origin = {-62, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
      Placement(visible = true, transformation(origin = {-34, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.room room(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {34, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
  	connect(prescribedHeatFlow.port, room.port_conv) annotation(
  		Line(points = {{-24, 30}, {0, 30}, {0, 32}, {24, 32}}, color = {191, 0, 0}));
  	connect(ramp.y, prescribedHeatFlow.Q_flow) annotation(
  		Line(points = {{-50, 30}, {-44, 30}}, color = {0, 0, 127}));
  	annotation(
      experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
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
  
  	A.pipe_ pipe_1(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {-102, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.pipe_ pipe_2(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {6, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.pipe_ pipe_3(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {6, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  	A.pipe_ pipe_4(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {-100, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  
  	Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 110, offset = 0, startTime = 5)  annotation(
  		Placement(visible = true, transformation(origin = {-66, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Medium, N_nominal = 1750, T_start = 273.15 + 30, checkValve = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.034, 0.04}, head_nominal = {39.0, 27.0, 22.8}), m_flow_start = 0.0001, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nParallel = 1, p_a_start = 101325, redeclare function powerCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticPower(V_flow_nominal = {0, 0.034, 0.04}, W_nominal = {5000, 14700, 17000}), use_N_in = true, use_T_start = true, use_powerCharacteristic = true)  annotation(
  		Placement(visible = true, transformation(origin = {-32, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0)  annotation(
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
  	replaceable package MediumAir = Media.Air.DryAirNasa;
  	inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.roomCycle roomCycle(redeclare package Medium = Medium, Water_initT = 273.15 + 30) annotation(
      Placement(visible = true, transformation(origin = {-28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 273.15 + 80, nPorts = 1, p = 101325)  annotation(
      Placement(visible = true, transformation(origin = {-94, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, nPorts = 1, p = 101325)  annotation(
      Placement(visible = true, transformation(origin = {-96, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.room room(redeclare package Medium = MediumAir) annotation(
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
  	A.HEX_HP hex_hp(InitialTemp = Water_initT, redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-56, 62}, extent = {{-12, 14}, {12, 36}}, rotation = 90)));
  	A.pipe_ pipe_4(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {-42, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  	A.pipe_ pipe_1(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {-56, 88}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  	A.pipe_ pipe_2(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
  		Placement(visible = true, transformation(origin = {10, 104}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
  	Modelica.Fluid.Vessels.ClosedVolume BT(redeclare package Medium = Medium, T_start = BT_initT, V = BT_vol, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 4, p_start = BT_initP, use_HeatTransfer = false, use_T_start = true, use_portsData = false)  annotation(
      Placement(visible = true, transformation(origin = {50, 76}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  	A.pipe_ pipe_3(redeclare package Medium = Medium, pipe_initialTemp = Water_initT) annotation(
      Placement(visible = true, transformation(origin = {22, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  	Modelica.Fluid.Machines.PrescribedPump pump(redeclare package Medium = Medium, N_nominal = 1000, T_start = Water_initT, V = 1, allowFlowReversal = false, checkValve = true, redeclare function efficiencyCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.constantEfficiency(eta_nominal = 0.9), energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, redeclare function flowCharacteristic = Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow(V_flow_nominal = {0, 0.25, 0.5}, head_nominal = {100, 60, 0}), m_flow_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, nParallel = 1, use_HeatTransfer = false, use_N_in = true, use_T_start = true)  annotation(
      Placement(visible = true, transformation(origin = {-12, -12}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  	Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 38.5, offset = 0, startTime = 5)  annotation(
      Placement(visible = true, transformation(origin = {-38, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Fluid.Sources.Boundary_ph boundary_beg(redeclare package Medium = Medium, nPorts = 1, p = 101325, use_h_in = true)  annotation(
  		Placement(visible = true, transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  	Modelica.Fluid.Sources.Boundary_pT boundary_end(redeclare package Medium = Medium, nPorts = 1, p = 101325)  annotation(
  		Placement(visible = true, transformation(origin = {-90, 34}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  	A.room room(redeclare package Medium = MediumAir) annotation(
  		Placement(visible = true, transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.roomCycle roomCycle(redeclare package Medium = Medium, Rad_initT = Water_initT, Water_initT = Water_initT) annotation(
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
  	replaceable package MediumAir = Media.Air.DryAirNasa;
  	inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	A.plant plant(redeclare package Medium = Medium, redeclare package MediumAir = MediumAir)  annotation(
  		Placement(visible = true, transformation(origin = {-2, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Blocks.Sources.Constant const(k = 273.15 + 20)  annotation(
  		Placement(visible = true, transformation(origin = {-84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Blocks.Math.Feedback feedback annotation(
  		Placement(visible = true, transformation(origin = {-56, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  	Modelica.Blocks.Continuous.PI pi(T = 300, initType = Modelica.Blocks.Types.Init.InitialState, k = 5)  annotation(
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
      Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
  end controller;
  annotation(
    uses(Modelica(version = "4.0.0"), Buildings(version = "9.1.0")));
end A;