% Given the amplitudes and phases of two cosine functions, this computes
% the phase at which their sum is is maximum.  The trig identity is
%
% a*sin(x+theta_a) + b*sin(x+theta_b) = c*sin(x+phi)                                        (1)
%
% where
%
% c^2 = a^2 + b^2 + 2*a*b*cos(theta_a-theta_b)                                              (2)
%
% tan(phi) = [a*sin(theta_a) + b*sin(theta_b)] / [a*cos(theta_a) + b*cos(theta_b)]          (3)
%
% I want (1) in terms of cosines. So use
%
% sin(z) = cos(z-pi/2)                                                                      (4)
%
%  
% Then (1) is
%
% a*sin(x+theta_a) + b*sin(x+theta_b) = c*sin(x+phi) =
%
% a*cos(x+theta_a-pi/2) + b*cos(x+theta_b-pi/2) =  c*sin(x+phi)                             (5)
%
% In my application, I only care about relative phases in the cosines, so set
% 
% theta_a-pi/2=0
%
% theta_a=pi/2                                                                              (6)
%
% I also only care about relative amplitudes of the cosines, so set 
%
% a=1                                                                                       (7)
%
% With (6,7),, Eq. (2) becomes
%
% c^2 = 1 + b^2 + 2*b*cos(pi/2 - theta_b)                                                   (8)
%
% and (3) becomes
%
% tan(phi) = [sin(pi/2) + b*sin(theta_b)] / [cos(pi/2) + b*cos(theta_b)]
%          = [1         + b*sin(theta_b)] /  b*cos(theta_b)                                 (9)
%
%
%
% and (5) becomes (a=1 & theta_a-pi/2=0)
%
% cos(x) + b*cos(x+theta_b-pi/2) =  c*sin(x+phi)                                           (10)
%
% Use (4) on right side to get
%
% cos(x) + b*cos(x+theta_b-pi/2) =  c*cos(x+phi-pi/2)                                      (11)
%
% 
% This program computes ?
% 
%
% INPUTS
%
% b                   1 x 1. See (5). For cases of usual interest, 0<=b<=1.
%
% OUTPUTS
%
% x                   1 x m.  Different values of x (rad)
% phi                 1 x m.  phi(i) yields maximum of c*cos(x+phi-pi/2).
% d_phase             1 x m.  Different values of Eq. (13)
% 
function []=sum_two_arbitrary_cosines(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progname='sum_two_arbitrary_cosines.m';

