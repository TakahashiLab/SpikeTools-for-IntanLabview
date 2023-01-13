%dn:data year
%an:animal label
function [homePath,dataPath]=PDdataC(dn,an)


%directory path,trigger band,stim LED,LFP left/right,Ctx tetrodes,Str tetrodes,led intensity,CC-INT,CC-PYR,Tag-PV(pos),Tag-Pyr(neg)
switch (dn)
  case '2018',
    homePath='/mnt/TakahashiLabA/mouse/2018/';
    %%%TK18052901
    switch (an)

        %%%%%
        %%%opto fail?
        %%%%%
      case 'TK18052901',%%PD%%OK
        dataPath={
            'TK18052901/2018061201','con',7,[1 6],[1:8],[],5,[13 15 29],[14 17 19 27],[],[]
            'TK18052901/2018061202','con',7,[1 6],[1:8],[],4,[],[],[],[]
            'TK18052901/2018061203','con',7,[1 6],[1:8],[],3,[],[],[],[]
            %            'TK18052901/2018061901','con',2,[1 6],[1:8],[],5,[18 19 21 28 33],[15 17 23 34 37 39 40],[],[]%%opto fail?
            %'TK18052901/2018061902','con',2,[1 6],[1:8],[],4,[],[],[],[]
            %            'TK18052901/2018061903','con',2,[1 6],[1:8],[],3,[],[],[],[]
                 };

      case 'TK18052902',%%PD%%OK
        dataPath={
            'TK18052902/2018060401','con',6,[1 5],[1:8],[],5,[ ],[],[],[]
            'TK18052902/2018060402','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18052902/2018060403','con',6,[1 5],[1:8],[],3,[],[],[],[]
            'TK18052902/2018060501','con',3,[1 5],[1:8],[],5,[4 15 ],[1 13 ],[],[]
            'TK18052902/2018060502','con',3,[1 5],[1:8],[],4,[],[],[],[] 
            'TK18052902/2018060503','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK18052902/2018061101','con',6,[1 5],[1:8],[],5,[6 ],[4 ],[24],[19]%%%clear tagging PV 24
            'TK18052902/2018061102','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18052902/2018061103','con',6,[1 5],[1:8],[],3,[],[],[],[]
                 };

      case 'TK18110801',%%PD control no stim%%OK
        dataPath={
            'TK18110801/2018112601','con',3,[1 5],[1:8],[],5,[3 4 7 11 26 32 45],[1 5 9 28 31 37 43],[],[]
            'TK18110801/2018112602','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK18110801/2018112603','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK18110801/2018112901','con',6,[1 5],[1:8],[],5,[4 8 16 17 29 32 33 36 38],[1 14 15 18 26 27 34 35 37],[38],[]
            'TK18110801/2018112902','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18110801/2018112903','con',6,[1 5],[1:8],[],3,[],[],[],[]
                 };

      case 'TK18121301',%%normal control%%OK
        dataPath={
            'TK18121301/2019010201','con',6,[1 5],[1:8],[],5,[9 37 50 60 61],[7 28 45 54 55 57],[],[]
            'TK18121301/2019010202','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18121301/2019010203','con',6,[1 5],[1:8],[],3,[],[],[],[]
            'TK18121301/2019010801','con',3,[1 5],[1:8],[],5,[7 10 41 54 55 60 63],[13 36 37 57 62],[],[8]
            'TK18121301/2019010802','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK18121301/2019010803','con',3,[1 5],[1:8],[],3,[],[],[],[]
                 };

      case 'TK18121401',%%PD%%OK
        dataPath={
            'TK18121401/2018121901','con',3,[1 5],[1:8],[],5,[2 ],[4],[],[]
            'TK18121401/2018121902','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK18121401/2018121903','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK18121401/2018123001','con',6,[1 5],[1:8],[],5,[21 31 34],[11 30 35],[],[39]
            'TK18121401/2018123002','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18121401/2018123003','con',6,[1 5],[1:8],[],3,[],[],[],[]
            'TK18121401/2019010701','con',6,[1 5],[1:8],[],5,[11 12 19 24 26 53 ],[1 13 15 16 17 18 22 50 52],[],[]
            'TK18121401/2019010702','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK18121401/2019010703','con',6,[1 5],[1:8],[],3,[],[],[],[]
                 };


      case 'TK19011801',%%PD%%OK
        dataPath={
            'TK19011801/2019020401','con',2,[1 5],[1:8],[],5,[],[],[],[]%checked
            'TK19011801/2019020402','con',2,[1 5],[1:8],[],4,[],[],[],[]
            'TK19011801/2019020403','con',2,[1 5],[1:8],[],3,[],[],[],[]
            'TK19011801/2019020701','con',7,[1 5],[1:8],[],5,[19],[],[19],[17]
            'TK19011801/2019020702','con',7,[1 5],[1:8],[],4,[],[],[],[]
            'TK19011801/2019020703','con',7,[1 5],[1:8],[],3,[],[],[],[]
                 };

      case 'TK19012501',%%PD%%OK
        dataPath={
            'TK19012501/2019021201','con',7,[1 5],[1:8],[],5,[5],[1],[],[17 29]%%nega-posi(inter-hemisphere) 5 inactive 1 active
            'TK19012501/2019021202','con',7,[1 5],[1:8],[],4,[],[],[],[]
            'TK19012501/2019021203','con',7,[1 5],[1:8],[],3,[],[],[],[]
            'TK19012501/2019021301','con',2,[1 5],[1:8],[],5,[14 10],[8 ],[],[]
            'TK19012501/2019021302','con',2,[1 5],[1:8],[],4,[],[],[],[]
            'TK19012501/2019021303','con',2,[1 5],[1:8],[],3,[],[],[],[]
                 };
    end
        %%%%%    20220810    
  case '2019',
    homePath='/mnt/TakahashiLabA/mouse/2019/';

    switch (an)
      case 'TK19040101',%%PD %%OK
        dataPath={
         'TK19040101/2019041001','con',3,[1 5],[1:8],[],5,[7 24 25 45 47],[1 26 27 28 46 47 48 49 50],[25 29 ],[21 23 27]%%25:PV-tagging BEST, 29:PV-tagging good,21,23: inhibition  BEST, contra: 45,47(pv):(initial inihibition?), 49(pyr): excitation?
            'TK19040101/2019041001','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19040101/2019041001','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19040101/2019041201','con',6,[1 5],[1:8],[],5,[ 7 11 16 28 46 49 51 52 55 57],[6 22 47 48 50 53],[38],[45 62]%%38:PV-tag Good
            'TK19040101/2019041201','con',6,[1 5],[1:8],[],3,[],[],[],[]
            'TK19040101/2019041707','con',6,[1 5],[1:8],[],5,[5 10 26 44 48 50],[1 11 12 25 42 45 46 47 52],[],[]%%44:rebound,26:contra-PV inhibition
            'TK19040101/2019041708','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19040101/2019041807','con',3,[1 5],[1:8],[],5,[3 18 40],[4 16 36 37 38 39],[5],[9 16 17 ]%%18 rebound-inhibition?, 35::contra-?? inhibition
            'TK19040101/2019041808','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19040101/2019041907','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19040101/2019041908','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19040101/2019042201','con',6,[1 5],[1:8],[],5,[16 17 29 32 33],[13 20 25 30 31],[28],[ 30]%%contra: 16(pv):inhibition, 20,25(pyr):(inhibition), ipsi: ,28(pv):excitation,30(pyr):inhibition
            'TK19040101/2019042202','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19040101/2019042203','con',6,[1 5],[1:8],[],3,[],[],[],[]
                 };

      case 'TK19050801',%%PD %%OK
        dataPath={
        %'TK19050801/2019051401','con',3,[1 5],[1:8],[],5,[20 23 30 31],[16 17 25 26 27 29],[18 20 21 23 31],[17 22 24 25 26 29]%%20:PV-tagging BEST2, 18:PV-tagging good, 21:PV mild
            'TK19050801/2019051402','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19050801/2019051403','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19050801/2019051507','con',6,[1 5],[1:8],[],5,[3 4 ],[2 7 8  ],[],[30]%%
            'TK19050801/2019051508','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19050801/2019051707','con',3,[1 5],[1:8],[],5,[9 11 36 39 ],[12 14 35],[],[]%%
            'TK19050801/2019051708','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19050801/2019051709','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19050801/2019052101','con',6,[1 5],[1:8],[],5,[],[],[],[]%%
            'TK19050801/2019052102','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19050801/2019052103','con',6,[1 5],[1:8],[],3,[],[],[],[]
                 };

        %%%%%
        %%%opto fail?
        %%%%%
      case 'TK19051601',%%PD%%OK
        dataPath={
            'TK19051601/2019052001','con',3,[1 5],[1:8],[],5,[32 33 34 35 36 ],[30 31 37],[],[]%%
            'TK19051601/2019052002','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19051601/2019052003','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19051601/2019052207','con',3,[1 5],[1:8],[],5,[8 24 30 56],[17 21 29 59 61],[],[]%%
            'TK19051601/2019052208','con',3,[1 5],[1:8],[],4,[],[],[],[]
            'TK19051601/2019052209','con',3,[1 5],[1:8],[],3,[],[],[],[]
            'TK19051601/2019052307','con',6,[1 5],[1:8],[],5,[7 21 34 38 39 61 63],[9 32 33 35 36 37 64],[],[]%%
            'TK19051601/2019052308','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19051601/2019052407','con',3,[1 5],[1:8],[],5,[7 10 11 17 23 45 58],[6 8 12 13 19 21 32 47 56],[],[]%%
            'TK19051601/2019052408','con',3,[1 5],[1:8],[],4,[],[],[],[]
                 };

      case 'TK19060301',%%PD 20220628%%%OK
        dataPath={
            'TK19060301/2019061101','con',6,[1 5],[1:8],[],5,[8 20 29 33 40],[4 18 27 34 38],[16 17],[]
            'TK19060301/2019061102','con',6,[1 5],[1:8],[],4,[],[],[],[]
            'TK19060301/2019061103','con',6,[1 5],[1:8],[],3,[],[],[],[]
            'TK19060301/2019061201','con',3,[1 5],[1:8],[],5,[1 4 7 8],[2 3],[1 4 7 11],[]%%%10:theta cell
            'TK19060301/2019061202','con',3,[1 5],[1:8],[],4,[],[],[],[]
                 };

    end
    
    case '2019f',
        homePath='f:/2019/';

        switch (an)
            
           case 'TK19081901',%PD, no ChRs on Str %%OK
                dataPath={
                    'TK19081901/2019100401','con',6,[1 5],[1:5 8],[6 7],5,[33],[31],[],[26]
                    'TK19081901/2019100402','con',6,[1 5],[1:5 8],[6 7],4,[],[],[],[]
                    'TK19081901/2019100403','con',6,[1 5],[1:5 8],[6 7],3,[],[],[],[]
                    'TK19081901/2019100404','con',3,[1 5],[1:5 8],[6 7],5,[10 22],[13 17],[5],[9 12 23]%%5:PV-tagging Best3
                    'TK19081901/2019100405','con',3,[1 5],[1:5 8],[6 7],4,[],[],[],[]
                    'TK19081901/2019100406','con',3,[1 5],[1:5 8],[6 7],3,[],[],[],[]
                         };

                %%very noisy
          case 'TK19100801',%PD,Str-Str(1 8), Ctx-LED %%OK
                dataPath={
                    'TK19100801/2019102301','con',7,[1 5],[3 4 5 6 7],[1 2 8],5,[],[],[],[]
                    'TK19100801/2019102302','con',7,[1 5],[3 4 5 6 7],[1 2 8],4,[],[],[],[]
                    'TK19100801/2019102303','con',7,[1 5],[3 4 5 6 7],[1 2 8],3,[],[],[],[]
                    'TK19100801/2019102304','con',2,[1 5],[3 4 5 6 7],[1 2 8],5,[],[],[],[]%%19: INT?
                    'TK19100801/2019102305','con',2,[1 5],[3 4 5 6 7],[1 2 8],4,[],[],[],[]
                    'TK19100801/2019102306','con',2,[1 5],[3 4 5 6 7],[1 2 8],3,[],[],[],[]
                    'TK19100801/2019110601','con',7,[1 5],[3 4 5 6 7],[1 2 8],5,[],[],[],[]
                    'TK19100801/2019110602','con',7,[1 5],[3 4 5 6 7],[1 2 8],4,[],[],[],[]
                    'TK19100801/2019110603','con',7,[1 5],[3 4 5 6 7],[1 2 8],3,[],[],[],[]
                    'TK19100801/2019110604','con',2,[1 5],[3 4 5 6 7],[1 2 8],5,[],[],[],[5 6]
                    'TK19100801/2019110605','con',2,[1 5],[3 4 5 6 7],[1 2 8],4,[],[],[],[]
                    'TK19100801/2019110606','con',2,[1 5],[3 4 5 6 7],[1 2 8],3,[],[],[],[]
                         };

          case 'AZ19091201',%normal %%OK no INT
                dataPath={
                    'AZ19091201/2019102101','con',6,[1 5],[1:8],[],5,[],[],[],[]
                    'AZ19091201/2019102102','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ19091201/2019102103','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ19091201/2019102104','con',3,[1 5],[1:8],[],5,[],[],[],[]%%9:inhibition rebound pyr?
                    'AZ19091201/2019102105','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ19091201/2019102106','con',3,[1 5],[1:8],[],3,[],[],[],[]
                         };

          case 'AZ19120501',%PD%%OK
                dataPath={
                    'AZ19120501/2020020601','con',3,[1 5],[1:8],[],5,[9 14],[6 11 12 13 15 ],[],[]
                    'AZ19120501/2020020602','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ19120501/2020020603','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ19120501/2020020604','con',6,[1 5],[1:8],[],5,[4 7 9],[6 8],[],[]%%12,rebound
                    'AZ19120501/2020020605','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ19120501/2020020606','con',6,[1 5],[1:8],[],3,[],[],[],[]
                         };

          case 'AZ20013101',%PD%%OK
                dataPath={
                    'AZ20013101/2020021901','con',6,[1 5],[1:8],[],5,[5 19 21],[4 17 18],[],[]%%no tagging
                    'AZ20013101/2020021902','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20013101/2020021903','con',6,[1 5],[1:8],[],3,[],[],[],[]
                         };

          case 'AZ20031201',%PD%% noisy
                dataPath={
                    'AZ20031201/2020032701','con',3,[1 5],[1:8],[],5,[13],[12],[],[]
                    'AZ20031201/2020032702','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031201/2020032703','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20031201/2020032704','con',6,[1 5],[1:8],[],5,[],[],[],[]
                    'AZ20031201/2020032705','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031201/2020032706','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20031201/2020040201','con',3,[1 5],[1:8],[],5,[18 22],[19 21],[9],[]
                    'AZ20031201/2020040202','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031201/2020040203','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20031201/2020040204','con',6,[1 5],[1:8],[],5,[3 ],[7],[],[]%4 contra-rebound, 18:contra-inhibittion
                    'AZ20031201/2020040205','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031201/2020040206','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    };

        end
    
  case '2020h',
        homePath='h:/2020/';

        switch (an)        
          case 'AZ20030601',%PD%%noisy
                dataPath={
                    'AZ20030601/2020040801','con',3,[1 5],[1:8],[],5,[],[],[],[]
                    'AZ20030601/2020040802','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20030601/2020040803','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20030601/2020040804','con',6,[1 5],[1:8],[],5,[4 19],[12 21],[13 22],[]%%13: good PV; 4, 6,12 contra-inhibition 
                    'AZ20030601/2020040805','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20030601/2020040806','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    };

          case 'AZ20031901',%PD%%OK
                dataPath={
                    'AZ20031901/2020041301','con',3,[1 5],[1:8],[],5,[],[],[],[4]%%
                    'AZ20031901/2020041302','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031901/2020041303','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20031901/2020041304','con',6,[1 5],[1:8],[],5,[5 30],[1 15 31 ],[37 40],[38 39 48]%%2,3,14: contra-excitation, 39:rebound-excitation
                    'AZ20031901/2020041305','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20031901/2020041306','con',6,[1 5],[1:8],[],3,[],[],[],[]
                         };
                
          case 'AZ20070101',%PD%%no spike
                dataPath={
                %                    'AZ20070101/2020072101','con',3,[1 5],[1:3 7 8],[4:6],5,[],[],[],[]%%nospike
                    'AZ20070101/2020072102','con',3,[1 5],[1:3 7 8],[4:6],4,[],[],[],[]
                    'AZ20070101/2020072103','con',3,[1 5],[1:3 7 8],[4:6],3,[],[],[],[]
                    'AZ20070101/2020072104','con',6,[1 5],[1:3 7 8],[4:6],5,[2],[1],[],[]%2:contra-excitation, 3:contra-inhibition
                    'AZ20070101/2020072105','con',6,[1 5],[1:3 7 8],[4:6],4,[],[],[],[]
                    'AZ20070101/2020072106','con',6,[1 5],[1:3 7 8],[4:6],3,[],[],[],[]
                    %                    'AZ20070101/2020072901','con',3,[1 5],[1:3 7 8],[4:6],5,[],[],[],[]%very noisy
                    'AZ20070101/2020072902','con',3,[1 5],[1:3 7 8],[4:6],4,[],[],[],[]
                    'AZ20070101/2020072903','con',3,[1 5],[1:3 7 8],[4:6],3,[],[],[],[]
                    'AZ20070101/2020072904','con',6,[1 5],[1:3 7 8],[4:6],5,[5],[8],[],[8]%%8:rebound
                    'AZ20070101/2020072905','con',6,[1 5],[1:3 7 8],[4:6],4,[],[],[],[]
                    'AZ20070101/2020072906','con',6,[1 5],[1:3 7 8],[4:6],3,[],[],[],[]
                    };

        end

  case '2020e',
        homePath='e:/2020/';

        switch (an)        
          case 'AZ20121501',%PD%%noisy
                dataPath={
                    'AZ20121501/2021011401','con',3,[1 5],[1:8],[],5,[],[],[],[]%%13: contra-inhibition
                    'AZ20121501/2021011402','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20121501/2021011403','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    %'AZ20121501/2021011404','con',6,[1 5],[1:8],[],5,[],[],[],[]%%no spike
                    'AZ20121501/2021011405','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20121501/2021011406','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20121501/2021012001','con',3,[1 5],[1:8],[],5,[6 11],[7 8 9],[],[]%%
                    'AZ20121501/2021012002','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20121501/2021012003','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    %'AZ20121501/2021012004','con',6,[1 5],[1:8],[],5,[],[],[],[]%%no spike
                    'AZ20121501/2021012006','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ20121501/2021012701','con',6,[1 5],[1:8],[],5,[5 13 16 19],[4 12 14],[],[]%%
                    'AZ20121501/2021012702','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ20121501/2021012703','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    };

          case 'AZ21012601',%PD%%OK
                dataPath={
                    'AZ21012601/2021020301','con',3,[1 5],[1:5 7:8],[6],5,[13 17 18 22],[9 15 16 20],[],[]%18 contra-inhibition
                    'AZ21012601/2021020302','con',3,[1 5],[1:5 7:8],[6],4,[],[],[],[]
                    'AZ21012601/2021020303','con',3,[1 5],[1:5 7:8],[6],3,[],[],[],[]
                    'AZ21012601/2021020304','con',6,[1 5],[1:5 7:8],[6],5,[12 13 14 21 34 35 38 41],[5 6 13 18 20 33 36 39 40 42 43],[35],[32 33 36]%%
                    'AZ21012601/2021020305','con',6,[1 5],[1:5 7:8],[6],4,[],[],[],[]
                    'AZ21012601/2021020306','con',6,[1 5],[1:5 7:8],[6],3,[],[],[],[]
                         };

          case 'AZ21010701',%PD%%poor
                dataPath={
                    'AZ21010701/2021030101','con',3,[1 5],[1:8],[],5,[5],[3],[],[]%5:contra-inhibition, 6:contra-rebound
                    'AZ21010701/2021030103','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ21010701/2021030106','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ21010701/2021030501','con',3,[1 5],[1:8],[],5,[22],[19],[],[]%15:contra-excitation
                    'AZ21010701/2021030502','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ21010701/2021030503','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ21010701/2021030504','con',6,[1 5],[1:8],[],5,[],[],[],[]%no information
                    'AZ21010701/2021030505','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ21010701/2021030506','con',6,[1 5],[1:8],[],3,[],[],[],[]
                         };

          case 'AZ21012501',%PD%%OK
                dataPath={
                    'AZ21012501/2021030801','con',3,[1 5],[1:2 4 5 7],[3 6 8],5,[8],[7 10],[],[7 10 11]
                    'AZ21012501/2021030802','con',3,[1 5],[1:2 4 5 7],[3 6 8],4,[],[],[],[]
                    'AZ21012501/2021030803','con',3,[1 5],[1:2 4 5 7],[3 6 8],3,[],[],[],[]

                    'AZ21012501/2021030804','con',6,[1 5],[1:2 4 5 7],[3 6 8],5,[1 2 5 6 11 15 17 25],[3 5 6 8 13 14 16 24],[],[28]%7,10,12,15,17,18,19,21: contra-rebound-excitation (INT?)
                    'AZ21012501/2021030805','con',6,[1 5],[1:2 4 5 7],[3 6 8],4,[],[],[],[]
                    'AZ21012501/2021030806','con',6,[1 5],[1:2 4 5 7],[3 6 8],3,[],[],[],[]
                         };
        end


  case '2022h',
        homePath='h:/2022/';

        switch (an)        
          case 'AZ22070601',%normal control without opt
                dataPath={
                    'AZ22070601/2022072501','con',3,[1 5],[1:8],[],5,[8 12 19],[5 11],[],[]
                    'AZ22070601/2022072502','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ22070601/2022072503','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ22070601/2022072504','con',6,[1 5],[1:8],[],5,[10 13 15],[6 7 14 17],[],[]
                    'AZ22070601/2022072505','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ22070601/2022072506','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    };
          case 'AZ22070701',%normal control without opt
                dataPath={
                    'AZ22070701/2022072801','con',3,[1 5],[1:8],[],5,[4 5 8 30 36 37],[2 7 24 38 39],[],[]%%no tagging OK?
                    'AZ22070701/2022072802','con',3,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ22070701/2022072803','con',3,[1 5],[1:8],[],3,[],[],[],[]
                    'AZ22070701/2022072804','con',6,[1 5],[1:8],[],5,[1 6 8 22],[2 5 7 9 10 13 23],[],[]
                    'AZ22070701/2022072805','con',6,[1 5],[1:8],[],4,[],[],[],[]
                    'AZ22070701/2022072806','con',6,[1 5],[1:8],[],3,[],[],[],[]
                    };
        end
end


