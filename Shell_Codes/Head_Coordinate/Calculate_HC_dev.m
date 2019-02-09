%% To calculate the subjects mean head position and their head coordinate deviation for each run (compared to their mean position)

% Create a index to store the fields and values (to create the final strucutre)
clear all
subj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
n = 1

for i = subj
    subject{n} = genvarname(['subject' int2str(i)])
    n = n + 1
    subject{n} = struct('mean',[],'S1R1',[],'S1R2',[],'S2R1',[], 'S2R2',[]);
    n = n + 1
end


%Load the head coordinate of each runs
n= 2
for i = subj
    eval(['subject' int2str(i)])
    head_filename = [subjectdata.subjectdir filesep 'head_coord' int2str(i) '.mat']
    load(head_filename)
    for S = [1 2]  %Pour les 2 sessions
        for R = [1 2]  %Pour les 2 runs
            if S == 1 &  R == 1
                name1 = head_coord.session1.run01
            elseif S == 1 &  R == 2
                name2 = head_coord.session1.run02
            elseif S == 2 & R == 1
                name3 = head_coord.session2.run01
            elseif S == 2 & R == 2
                name4 = head_coord.session2.run02
            end
        end
    end
    
    % Calculate the mean position
    for colomn = 1:3
    for ligne = 1:3
    X= [name1(ligne,colomn); name2(ligne,colomn); name3(ligne,colomn); name4(ligne,colomn)]
    subject{1, n}.mean(ligne,colomn) = mean(X, 1)
    end
    end
    
     % Calculate the deviation for each run
    k = 1 
    session = 1 
    run = 1
    
    for runs = 1:4
    for colomn = 1 : 3
    for ligne = 1:3
    name{k} = genvarname(['name' int2str(k)])
    name{k} = eval(name{k})
            if name{k} == name1
                subject{1, n}.S1R1(ligne,colomn) = abs(name{k}(ligne,colomn) - subject{1, n}.mean(ligne,colomn))
            elseif name{k} == name2
                subject{1, n}.S1R2(ligne,colomn) = abs(name{k}(ligne,colomn) - subject{1, n}.mean(ligne,colomn))
            elseif name{k} == name3
                subject{1, n}.S2R1(ligne,colomn) = abs(name{k}(ligne,colomn) - subject{1, n}.mean(ligne,colomn))
            elseif name{k} == name4
                subject{1, n}.S2R2(ligne,colomn) = abs(name{k}(ligne,colomn) - subject{1, n}.mean(ligne,colomn))
            end
    end
    end
    k = k +1
    end
    n = n + 2    
    clear('head_coord')
end

% Store the calulations for each subjects in a structure 
head_coord_dev = struct(subject{:})   

% Display the results
subj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
for i = subj
 disp(['head_coord_dev.subject' int2str(i) '.S1R1']) 
 eval(sprintf(['head_coord_dev.subject' int2str(i) '.S1R1'])) 
 disp(['head_coord_dev.subject' int2str(i) '.S1R2']) 
 eval(sprintf(['head_coord_dev.subject' int2str(i) '.S1R2'])) 
 disp(['head_coord_dev.subject' int2str(i) '.S2R1']) 
 eval(sprintf(['head_coord_dev.subject' int2str(i) '.S2R1'])) 
 disp(['head_coord_dev.subject' int2str(i) '.S2R2']) 
 eval(sprintf(['head_coord_dev.subject' int2str(i) '.S2R2'])) 
end




