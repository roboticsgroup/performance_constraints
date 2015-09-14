function lwr = lwr()

%Creates a KUKA LWR IV for Robotics Toolbox with DH parameters. 
%Returns robot struct.
L1 = Link('d', 0.3105, 'a', 0, 'alpha', pi/2, 'offset', 0);
L2 = Link('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0);
L3 = Link('d', 0.4, 'a', 0, 'alpha',-pi/2, 'offset', 0);
L4 = Link('d', 0, 'a', 0, 'alpha', pi/2, 'offset', 0);
L5 = Link('d', 0.39, 'a', 0, 'alpha', pi/2);
L6 = Link('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0);
L7 = Link('d', 0.078, 'a', 0, 'alpha', 0, 'offset', 0);
lwr = SerialLink([L1 L2 L3 L4 L5 L6 L7], 'name', 'LWR 7DOF');

end