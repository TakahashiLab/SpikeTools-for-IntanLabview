function batchBhv(opt)
if nargin==0
    opt='win';
end

if opt=='win'
    yn='2018';
    dn={
        'TK18052901',
        'TK18052902',%%PD    
        'TK18110801',%%PD control no stim
        'TK18121301',%%normal control
        'TK18121401',%%PD
        'TK19011801',%%PD
        'TK19012501',%%PD
       };


    for i=1:size(dn,1)
        bhv=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'bhv');
    end



    yn='2019';
    dn={'TK19040101',%%PD 
        'TK19050801',%%PD 
        'TK19051601',%%PD
        'TK19060301',%%PD 20220628
       };

    for i=1:size(dn,1)
        bhv=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'bhv');
    end
else
    
yn='2019f';

dn={
    'TK19081901',%PD, no ChRs on Str
    'TK19100801',%PD,Str-Str(1 8), Ctx-LED
    'AZ19091201',%normal
    'AZ19120501',%PD
    'AZ20013101',%PD
    'AZ20031201',%PD
};

    for i=1:size(dn,1)
        bhv=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'bhv');
    end
    
yn='2020h';

dn={
    'AZ20030601',%PD
    'AZ20031901',%PD
     'AZ20070101'%PD
    };

    for i=1:size(dn,1)
        bhv=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'bhv');
    end


yn='2020e';
dn={'AZ20121501',%PD
     'AZ21012601',%PD
     'AZ21010701',%PD
     'AZ21012501'%PD
};

    for i=1:size(dn,1)
        bhv=proPD(yn,dn{i},'method','behavior');
        save([dn{i} '.mat'],'bhv');
    end
end
return;


