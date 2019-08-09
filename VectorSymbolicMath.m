%Written By: Dean Keithly
%Written On: 7/29/2019

%% Define Symbols
syms r wpv ii a e W w E v mu X Y Z Xdot Ydot Zdot


%% Define Vectors
% reqn1 = r == (a*(1-e^2))/(1+e*(cos(E)-e)/(1-e*cos(E)))
% reqn2 = r == (X^2+Y^2+Z^2)^0.5
% 
% wpveqn = wpv == w+v
% 
% Xeqn = X == r*(cos(W)*cos(wpv)-sin(W)*sin(wpv)*cos(ii))
% Yeqn = Y == r*(sin(W)*cos(wpv)-cos(W)*sin(wpv)*cos(ii))
% Zeqn = Z == r*sin(ii)*sin(wpv)
% 
% rv = [X,Y,Z]  
%   
% rdotPQW = [-(mu/a)^0.5*(sin(E))/(1-e*cos(E))
%     (mu/(a*(1-e^2)))^0.5*(e+(cos(E)-e)/(1-e*cos(E)))
%     0]
% 
% %rotation matrix between the two
% ijkQpqw = [cos(W)*cos(w)-sin(W)*sin(w)*cos(ii), -cos(W)*sin(w)-sin(W)*cos(w)*cos(ii), -sin(W)*sin(ii);
%     sin(W)*cos(w)+cos(W)*sin(w)*cos(ii), -sin(W)*sin(w)+cos(W)*cos(w)*cos(ii), -cos(W)*sin(ii);
%     sin(w)*sin(ii), cos(w)*sin(ii), cos(ii)]
% 
% %Calculate rdot in ijk
% rdotijk = ijkQpqw*rdotPQW
% 
% Xdoteqn = Xdot == rdotijk(1)
% Ydoteqn = Ydot == rdotijk(2)
% Zdoteqn = Zdot == rdotijk(3)
% 
% h = cross(rv,rdotijk)
% ev = (cross(rdotijk,h))/mu-rv/norm(rv) 
% eeqn = e == norm(ev)
% 
% %n from https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
% n = [-h(2), h(1), 0]
% 
% ieqn = ii == acos(h(3))
% Weqn = W == acos(n(1)/norm(n)) % for ny > 0
% 
% weqn = w == acos(dot(n,ev)/(norm(ev)*norm(n)))
% 
% veqn = v == acos(dot(ev,rv)/(norm(ev)*norm(rv)))
% 
% %eqns = [reqn1,reqn2,weqn,Weqn,ieqn,eeqn, Xeqn, Yeqn, Zeqn, Xdoteqn, Ydoteqn, Zdoteqn]
% eqns = [reqn2, weqn, Weqn, ieqn, wpveqn, veqn, Xeqn, Yeqn, Zeqn, Xdoteqn, Ydoteqn, Zdoteqn]
% vars = [r,W,wpv,ii]
% 
% out = solve(eqns,vars)

%% From RV2COE in Vallado
%syms real X Y Z Xd Yd Zd p a e ii W w v mu
X = sym('X','real');
Y = sym('Y','real');
Z = sym('Z','real');
Xd = sym('Xd','real');
Yd = sym('Yd','real');
Zd = sym('Zd','real');
p = sym('p','real');
a = sym('a','real');
e = sym('e','real');
ii = sym('ii','real');
W = sym('W','real');
w = sym('w','real');
v = sym('v','real');
mu = sym('mu','real');


rvect = [X Y Z]
rscalar = norm(rvect)
rdvect = [Xd Yd Zd]
rdscalar = norm(rdvect)
hvect = cross(rvect,rdvect)
hscalar = norm(hvect)

nvect = cross([0 0 1],hvect)

evect = (((rdscalar^2-mu/rscalar)*rvect)-dot(rvect,rdvect)*rdvect)/mu

eeqn = e == norm(evect)
ieqn = ii == acos(hvect(3)/hscalar)
Weqn = W == acos(nvect(1)/norm(nvect))
weqn = w == acos(dot(nvect,evect)/norm(nvect)/norm(evect))
veqn = v == acos(dot(evect,rvect)/norm(evect)/norm(rvect))
aeqn = a == -(mu/2)/((rdscalar^2)/2-(mu/rscalar))
%peqn = p == a*(1-e^2)
reqn = rscalar == a*(1-e^2)/(1+e*cos(v))

eqns = [eeqn, ieqn, Weqn, weqn, veqn, aeqn, reqn]
vars = [a,e,ii,W,w,v]
out = solve(eqns,vars)

%To access representation, out.a, out.e, ...

%% Calculate Partial Derivatives
a11 = diff(rscalar,X)
a12 = diff(rscalar,Y)
a13 = diff(rscalar,Z)

a21 = diff(out.w+out.v,X)
a22 = diff(out.w+out.v,Y)
a23 = diff(out.w+out.v,Z)

a31 = diff(out.W,X)
a32 = diff(out.W,Y)
a33 = diff(out.W,Z)

a41 = diff(out.ii,X)
a42 = diff(out.ii,Y)
a43 = diff(out.ii,Z)


