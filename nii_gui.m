function nii_gui
%sample analysis of an event related design

prompt = {'Number of conditions:','Event duration (seconds, 0=variable):','Number of fMRI sessions:','TR (seconds):','Slice order (0=auto,1=[1234],2=[4321],3=[1324],4=[4231],5=[2143]:',...
'Use phasemap (1=yes)'};
dlg_title = 'Options for analyzing ';
def = {'2','1','1','0','0','0'};
answer = inputdlg(prompt,dlg_title,1,def);
if isempty(answer), return; end;
nCond = str2double(answer{1});
p.duration{1} = str2double(answer{2});
isVariableDuration = (p.duration{1,1} == 0);
nSess = str2double(answer{3});
p.TRsec = str2double(answer{4});
p.slice_order = str2double(answer{5});
isPhase = str2double(answer{6});
p.fmriname = [];
for s = 1: nSess
   p.fmriname = strvcat(p.fmriname, spm_select(1,'image',sprintf('Select 1st vol of session %d',s))); %#ok<REMFF1>
end
p.t1name = spm_select(1,'image','Select T1 image');
if isPhase
    p.phase =  spm_select(1,'image','Select fieldmap phase image');
    p.magn = spm_select(1,'image','Select fieldmap magnitude image');
else
    p.phase =  ''; %phase image from fieldmap
    p.magn = ''; %magnitude image from fieldmap    
end
for c = 1: nCond
    p.names{c} = ['Condition' num2str(c)];
    for s = 1: nSess
        dlg_title = ['Session ' num2str(s) ' Condition ' num2str(c)];
        prompt = {'Condition name', 'Onset times (sec):'};
        def = {p.names{c},''};
        if isVariableDuration
            prompt{3} = 'Event durations (sec)';
            def{3} = '';
            answer = inputdlg(prompt,dlg_title,[1 60],def);
            p.duration{s,c}= str2num(answer{3}); %#ok<ST2NM>
        else
            answer = inputdlg(prompt,dlg_title,[1 60],def);
        end
        p.names{c} = answer{1};
        p.onsets{s,c}= str2num(answer{2}); %#ok<ST2NM>
    end %for s: session
end %for c: condition
%statistical information (optional: if not provided not statitics)
[pth nam] = spm_fileparts(p.fmriname(1,:));
save(fullfile(pth, [nam '_batch']), '-struct', 'p');
%next: run the analysis...
nii_batch12(p);
