%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRAVEHEART - Open source software for electrocardiographic and vectorcardiographic analysis
% delete_beat_GUI.m -- Part of BRAVEHEART GUI
% Copyright 2016-2024 Hans F. Stabenau and Jonathan W. Waks
% 
% Source code/executables: https://github.com/BIVectors/BRAVEHEART
% Contact: braveheart.ecg@gmail.com
% 
% BRAVEHEART is free software: you can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
%
% BRAVEHEART is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. 
% If not, see <https://www.gnu.org/licenses/>.
%
% This software is for research purposes only and is not intended to diagnose or treat any disease.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delete_beat_GUI(hObject, eventdata, handles)

ecg = handles.ecg;
vcg = handles.vcg;
beats = handles.beats;
aps = pull_guiparams(hObject, eventdata, handles); 

set(handles.success_txt,'Visible','Off')

% Assign chosen beat as active_beat_number
str_index = get(handles.activebeats_list,'Value');
handles.active_beat_number = str_index;

str_matrix = get(handles.activebeats_list,'String');

% Delete beat from beats class
beats = beats.delete(str_index);
handles.beats = beats;  % update beats in handles
guidata(hObject, handles);

listbox_beats = beats_to_listbox(beats.Q, beats.QRS, beats.S, beats.Tend);
set(handles.activebeats_list,'String',listbox_beats);

    
% Finds outliers 
    beats = beats.find_outliers(vcg,aps);
    handles.beats = beats;
    guidata(hObject, handles);  % update handles
    
% Finds PVCs
    beats = beats.find_pvcs(vcg,aps);
    handles.beats = beats;
    guidata(hObject, handles);  % update handles


 
% Mark outliers and PVCs in listbox and mark in rhythm leads    
    graph_outliers_pvcs(beats, vcg, aps, hObject,eventdata,handles)
    
    
end