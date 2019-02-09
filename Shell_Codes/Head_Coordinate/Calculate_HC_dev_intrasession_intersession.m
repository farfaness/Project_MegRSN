%% To calculate the subjects mean head position for each run and their head coordinate deviation intra session (compared to their mean position) and inter session

% Create a index to store the fields and values (to create the final strucutre)
clear all
subj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
n = 1

for i = subj
    subject{n} = genvarname(['subject' int2str(i)])
    n = n + 1
    subject{n} = struct('mean_S1',[],'mean_S2',[],'mvt_intra_session1',[],'mvt_intra_session2',[],'mvt_inter_session',[]);
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
    
    % Calculate the mean position for each session
    for colomn = 1:3
    for ligne = 1:3
    session1= [name1(ligne,colomn); name2(ligne,colomn)]
    subject{1, n}.mean_S1(ligne,colomn) = mean(session1, 1)
    end
    end
    
    for colomn = 1:3
    for ligne = 1:3
    session2= [name3(ligne,colomn); name4(ligne,colomn)]
    subject{1, n}.mean_S2(ligne,colomn) = mean(session2, 1)
    end
    end
    
     % Calculate the mouvement intra_session1 for the first session
    subject{1, n}.mvt_intra_session1(1,1) = sqrt((name1(1,1) - name2(1,1))^2 + (name1(2,1) - name2(2,1))^2 + (name1(3,1) - name2(3,1))^2)
    subject{1, n}.mvt_intra_session1(1,2) = sqrt((name1(1,2) - name2(1,2))^2 + (name1(2,2) - name2(2,2))^2 + (name1(3,2) - name2(3,2))^2)
    subject{1, n}.mvt_intra_session1(1,3) = sqrt((name1(1,3) - name2(1,3))^2 + (name1(2,3) - name2(2,3))^2 + (name1(3,3) - name2(3,3))^2)
 
     % Calculate the mouvement intra_session2 for the second session
    subject{1, n}.mvt_intra_session2(1,1) = sqrt((name3(1,1) - name4(1,1))^2 + (name3(2,1) - name4(2,1))^2 + (name3(3,1) - name4(3,1))^2)
    subject{1, n}.mvt_intra_session2(1,2) = sqrt((name3(1,2) - name4(1,2))^2 + (name3(2,2) - name4(2,2))^2 + (name3(3,2) - name4(3,2))^2)
    subject{1, n}.mvt_intra_session2(1,3) = sqrt((name3(1,3) - name4(1,3))^2 + (name3(2,3) - name4(2,3))^2 + (name3(3,3) - name4(3,3))^2)
    
    % Calculate the mouvement inter_session
    
    subject{1, n}.mvt_inter_session(1,1) = sqrt((subject{1, n}.mean_S1(1,1) - subject{1, n}.mean_S2(1,1))^2 + (subject{1, n}.mean_S1(2,1) - subject{1, n}.mean_S2(2,1))^2 + (subject{1, n}.mean_S1(3,1) - subject{1, n}.mean_S2(3,1))^2)
    subject{1, n}.mvt_inter_session(1,2) = sqrt((subject{1, n}.mean_S1(1,2) - subject{1, n}.mean_S2(1,2))^2 + (subject{1, n}.mean_S1(2,2) - subject{1, n}.mean_S2(2,2))^2 + (subject{1, n}.mean_S1(3,2) - subject{1, n}.mean_S2(3,2))^2)
    subject{1, n}.mvt_inter_session(1,3) = sqrt((subject{1, n}.mean_S1(1,3) - subject{1, n}.mean_S2(1,3))^2 + (subject{1, n}.mean_S1(2,3) - subject{1, n}.mean_S2(2,3))^2 + (subject{1, n}.mean_S1(3,3) - subject{1, n}.mean_S2(3,3))^2)
  
    n = n + 2    
    clear('head_coord')
end


% Store the calulations in a structure 
head_coord_dev_intra_inter_session = struct(subject{:})   

% Display the results
subj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
for i = subj
 disp(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session1']) 
 eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session1' '(1,1)'])) 
 disp(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session2']) 
 eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session2'])) 
 disp(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_inter_session']) 
 eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_inter_session'])) 
end

% Save the results in a matrice

A = NaN(20, 9)
% A(1,1) = printmat('Subject')
% A(1,2) = 'mvt_intra_session1-leftear'
% A(1,3) = 'mvt_intra_session1-nasion'
% A(1,4) = 'mvt_intra_session1-rightear'
% A(1,5) = 'mvt_intra_session2-leftear'
% A(1,6) = 'mvt_intra_session2-nasion'
% A(1,7) = 'mvt_intra_session2-rightear'
% A(1,8) = 'mvt_inter_session-leftear'
% A(1,9) = 'mvt_inter_session-nasion'
% A(1,10) = 'mvt_inter_session-rightear'

subj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
n = 1

for i = subj
% A(n, 1) = i  
A(n,1) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session1' '(1,1)'])) 
A(n,2) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session1' '(1,2)'])) 
A(n,3) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session1' '(1,3)'])) 
A(n,4) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session2' '(1,1)'])) 
A(n,5) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session2' '(1,2)'])) 
A(n,6) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_intra_session2' '(1,3)'])) 
A(n,7) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_inter_session' '(1,1)'])) 
A(n,8) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_inter_session' '(1,2)'])) 
A(n,9) = eval(sprintf(['head_coord_dev_intra_inter_session.subject' int2str(i) '.mvt_inter_session' '(1,3)'])) 

n = n+1
end

printmat(A, 'Mvt', '1 2 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22', 'mvt_intra_session1-leftear mvt_intra_session1-nasion mvt_intra_session1-rightear mvt_intra_session2-leftear mvt_intra_session2-nasion mvt_intra_session2-rightear mvt_inter_session-leftear mvt_inter_session-nasion mvt_inter_session-rightear') 





