function [config, store] = nmfthgrInit(config)                            
% nmfthgrInit INITIALIZATION of the expLanes experiment NMFThresholdGrafic
%    [config, store] = nmfthgrInit(config)                                
%      - config : expLanes configuration state                            
%      -- store  : processing data to be saved for the other steps        
                                                                          
% Copyright: <userName>                                                   
% Date: 21-Dec-2017                                                       
                                                                          
if nargin==0, NMFThresholdGrafic(); return; else store=[];  end           
