%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by:                                          %
% - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
%                                                                           %
% Please send your comments to: amirzend@ku.edu                             %
%                                                                           %
% The code is intended for educational purposes and theoretical details     %
% are discussed in the textbook:                                            %
% Introduction to Shape and Topology Optimization using MATLAB              %
%                                                                           %
% Disclaimer:                                                               %
% The authors reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event        %
% caused by the use of the program.                                         %
%                                                                           %
% License:                                                                  %
% This software is used, copied and distributed under the licensing         %
% agreement contained in the file LICENSE in the top directory of           %
% the distribution.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath(cd))
format compact;close all; clear; clc
disp('STORX: Shape and Topology Optimization for Research and Experimentation');
disp('Version 2026.03');
disp('Copyright: Amir M. Mirzendehdel (amirzend@ku.edu) and Krishnan Suresh (ksuresh@wisc.edu)');
thisDir = cd;
%% Latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultBubblelegendInterpreter','latex');  
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultConstantlineInterpreter','latex');  
set(groot, 'defaultGraphplotInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter','latex');  
set(groot, 'defaultTextInterpreter','latex');  
set(groot, 'defaultTextarrowshapeInterpreter','latex');
set(groot, 'defaultTextboxshapeInterpreter','latex');

%% Font size
DEFAULT_FONT_SIZE = 18;
set(groot,'defaultAxesFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultBubblelegendFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultColorbarFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultConstantlineFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultGeoaxesFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultGraphplotEdgeFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultGraphplotNodeFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultLegendFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultPolaraxesFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultTextFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultTextarrowshapeFontSize',DEFAULT_FONT_SIZE)
set(groot,'defaultTextboxshapeFontSize',DEFAULT_FONT_SIZE)
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.0);
%% Warning
warning('off','all')