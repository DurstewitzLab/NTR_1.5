function array_vals=config_values(n_confs,min_value,max_value,is_integer)
%Auxiliar file for speeding-up parallel processing, returns a value in between min_value
%and max_value, equally or logaritmically spaced
%
%This code is not commercial and thus is not guaranteed but is used on regular basis in 
%multiple datasets.
%Based on: Balaguer-Ballester E, Lapish C, Seamans JK and Durstewitz D 2011
% "Attracting Dynamics of Frontal Cortex Populations during Memory Guided Decision-Making".
% PLoS Comput. Biol.
%Distributed accourding to the General Public LIcense GNU 3.0
%http://www.gnu.org/copyleft/gpl.html
%
%by Emili Balaguer-Ballester.
%Contact adresses: eb-ballester@bournemouth.ac.uk
%                  daniel.durstewitz@zi-mannheim.de
%                  emili.balaguer@uv.es
%
%
%Log of changes: Please specify here your code modifications
%(#Line@author: Change description). This will be useful if assistance is
%required
%
%
%__________________________________________________________________________
%
if nargin<6, is_integer=0; end; %if nargin<5, is_log=0; end;
%
step=(max_value-min_value)/n_confs;
if is_integer, step=floor(step); end 
if min_value+step>max_value, error('Maximum value is too small'), end
%
if step>0
array_vals=(min_value:step:max_value);
else
    array_vals=min_value;
end
if length(array_vals)<n_confs, 
    array_vals=[array_vals,step+array_vals(end)]; 
elseif length(array_vals)>n_confs, array_vals=array_vals(1:n_confs);
end


