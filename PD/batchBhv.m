function batchBhv(opt)
if nargin==0
    opt='win';
end

if opt=='win'
    yn='2018';
    dn={
        'TK18052901',
        'TK18052902',%%PD    
        'TK18121401',%%PD
        'TK19011801',%%PD
        'TK19012501',%%PD
       };


    for i=1:size(dn,1)
        [LED,L,R,~,details]=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'LED','L','R','details');
    end

    yn='2018';
    dn={
        'TK18110801',%%PD control no stim
        };
    [LED,L,R,~,details]=proPD(yn,dn{1},'method','behavior');
    save(['pdWOstim/' dn{1} '.mat'],'LED','L','R','details');


    yn='2018';
    dn={
        'TK18121301',%%normal control
        };
    [LED,L,R,~,details]=proPD(yn,dn{1},'method','behavior');
    save(['normal/' dn{1} '.mat'],'LED','L','R','details');

    yn='2019';
    dn={'TK19040101',%%PD 
        'TK19050801',%%PD 
        'TK19051601',%%PD
        'TK19060301',%%PD 20220628
       };

    for i=1:size(dn,1)
        [LED,L,R,~,details]=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'LED','L','R','details');
    end
else
    
yn='2019f';
dn={
    'TK19100801',%PD,Str-Str(1 8), Ctx-LED
    'AZ19120501',%PD
    'AZ20013101',%PD
    'AZ20031201',%PD
};

    for i=1:size(dn,1)
        [LED,L,R,~,details]=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'LED','L','R','details');    'TK19081901',%PD, no ChRs on Str

    end

yn='2019f';
dn={
    'AZ19091201',%normal
};
    [LED,L,R,~,details]=proPD(yn,dn{1},'method','behavior');
    save(['normal/' dn{1} '.mat'],'LED','L','R','details');

yn='2019f';
dn={
    'TK19081901',%PD, no ChRs on Str    
};
    [LED,L,R,~,details]=proPD(yn,dn{1},'method','behavior');
    save(['LEDquestion/' dn{1} '.mat'],'LED','L','R','details');

yn='2020h';

dn={
    'AZ20030601',%PD
    'AZ20031901',%PD
    };

    for i=1:size(dn,1)
        [LED,L,R,~,details]=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'LED','L','R','details');

    end

yn='2020h';

dn={
     'AZ20070101'%PD
};
    [LED,L,R,~,details]=proPD(yn,dn{1},'method','behavior');
    save(['ctx-str/' dn{1} '.mat'],'LED','L','R','details');

yn='2020e';
dn={'AZ20121501',%PD
     'AZ21012601',%PD
     'AZ21010701',%PD
     'AZ21012501'%PD
};

    for i=1:size(dn,1)
        [LED,L,R,~,details]=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'LED','L','R','details');
    end
end
return;





