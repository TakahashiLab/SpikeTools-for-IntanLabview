%dn:data year
%an:animal label
function [homePath,dataPath]=PDdata(dn,an)


%directory path,trigger band,stim LED,LFP,Ctx tetrodes,Str tetrodes
switch (dn)
  case '2018',
    homePath='/mnt/TakahashiLabA/mouse/2018/';
    %%%TK18052901
    switch (an)
      case 'TK18052901',%%PD
        dataPath={
            'TK18052901/2018061204','beta',7,1,[1:8],[]
            'TK18052901/2018061205','theta',7,1,[1:8],[]
            'TK18052901/2018061301','gamma',7,1,[1:8],[]
            'TK18052901/2018061302','gamma',7,8,[1:8],[]
            'TK18052901/2018061303','theta',7,8,[1:8],[]
            'TK18052901/2018061304','beta',7,8,[1:8],[]
            'TK18052901/2018061305','beta',2,8,[1:8],[]
            'TK18052901/2018061306','theta',2,8,[1:8],[]
            'TK18052901/2018061307','gamma',2,8,[1:8],[]
            'TK18052901/2018061308','beta',2,1,[1:8],[]
            'TK18052901/2018061501','beta',2,1,[1:8],[]
            'TK18052901/2018061502','theta',2,1,[1:8],[]
            'TK18052901/2018061503','gamma',2,1,[1:8],[]
            'TK18052901/2018061904','beta',7,3,[1:8],[]
            'TK18052901/2018061905','theta',7,3,[1:8],[]
            'TK18052901/2018061906','gamma',7,3,[1:8],[]
            'TK18052901/2018061907','beta',7,7,[1:8],[]
            'TK18052901/2018061908','theta',7,7,[1:8],[]
            'TK18052901/2018061909','gamma',7,7,[1:8],[]
            'TK18052901/2018062201','beta',2,5,[1:8],[]
            'TK18052901/2018062202','theta',2,5,[1:8],[]
            'TK18052901/2018062203','gamma',2,5,[1:8],[]
                 };

      case 'TK18052902',%%PD
        dataPath={
            'TK18052902/2018060404','beta',6,7,[1:8],[]
            'TK18052902/2018060405','theta',6,7,[1:8],[]
            'TK18052902/2018060406','gamma',6,7,[1:8],[]
            'TK18052902/2018060504','beta',3,7,[1:8],[]
            'TK18052902/2018060505','theta',3,7,[1:8],[]
            'TK18052902/2018060601','beta',6,1,[1:8],[]
            'TK18052902/2018060602','theta',6,1,[1:8],[]
            'TK18052902/2018060603','gamma',6,1,[1:8],[]
            'TK18052902/2018060604','beta',6,8,[1:8],[]
            'TK18052902/2018060605','theta',6,8,[1:8],[]
            'TK18052902/2018060606','gamma',6,8,[1:8],[]
            'TK18052902/2018060701','beta',3,4,[1:8],[]
            'TK18052902/2018060702','theta',3,4,[1:8],[]
            'TK18052902/2018060703','gamma',3,4,[1:8],[]
            'TK18052902/2018060704','beta',3,8,[1:8],[]
            'TK18052902/2018060705','theta',3,8,[1:8],[]
            'TK18052902/2018060706','gamma',3,8,[1:8],[]
            'TK18052902/2018062101','beta',6,1,[1:8],[]
            'TK18052902/2018062102','theta',6,1,[1:8],[]
            'TK18052902/2018062103','gamma',6,1,[1:8],[]
            'TK18052902/2018062104','beta',6,7,[1:8],[]
            'TK18052902/2018062105','theta',6,7,[1:8],[]
            'TK18052902/2018062106','gamma',6,7,[1:8],[]
                 };

      case 'TK18110801',%%PD control no stim
        dataPath={
            'TK18110801/2018112801','beta',3,4,[1:8],[]
            'TK18110801/2018112802','theta',3,4,[1:8],[]
            'TK18110801/2018112803','gamma',3,4,[1:8],[]
            'TK18110801/2018112804','beta',3,7,[1:8],[]
            'TK18110801/2018112805','theta',3,7,[1:8],[]
            'TK18110801/2018112806','gamma',3,7,[1:8],[]
            'TK18110801/2018112904','beta',6,5,[1:8],[]
            'TK18110801/2018112905','theta',6,5,[1:8],[]
            'TK18110801/2018112906','gamma',6,5,[1:8],[]
            'TK18110801/2018120301','beta',6,4,[1:8],[]
            'TK18110801/2018120302','theta',6,4,[1:8],[]
            'TK18110801/2018120303','gamma',6,4,[1:8],[]
                 };

case 'TK18121301',%%normal control,ChR2
        dataPath={
            'TK18121301/2019010401','beta',6,2,[1:8],[]
            'TK18121301/2019010402','theta',6,2,[1:8],[]
            'TK18121301/2019010403','gamma',6,2,[1:8],[]
            'TK18121301/2019010404','beta',6,7,[1:8],[]
            'TK18121301/2019010405','theta',6,7,[1:8],[]
            'TK18121301/2019010406','gamma',6,7,[1:8],[]
            'TK18121301/2019010701','beta',3,2,[1:8],[]
            'TK18121301/2019010702','theta',3,2,[1:8],[]
            'TK18121301/2019010703','gamma',3,2,[1:8],[]
            'TK18121301/2019010704','beta',3,7,[1:8],[]
            'TK18121301/2019010705','theta',3,7,[1:8],[]
            'TK18121301/2019010804','gamma',3,7,[1:8],[]
            'TK18121301/2019010805','theta',3,7,[1:8],[]
            %            'TK18121301/2019010806','beta',3,7,[1:8],[]
                 };

      case 'TK18121401',%%PD
        dataPath={
            'TK18121401/2018122801','beta',3,2,[1:8],[]
            'TK18121401/2018122802','theta',3,2,[1:8],[]
            'TK18121401/2018122803','gamma',3,2,[1:8],[]
            'TK18121401/2018122804','beta',3,7,[1:8],[]
            'TK18121401/2018122805','theta',3,7,[1:8],[]
            'TK18121401/2018122806','gamma',3,7,[1:8],[]
            'TK18121401/2018123004','beta',6,2,[1:8],[]
            %            'TK18121401/2018123005','theta',6,2,[1:8],[]
            'TK18121401/2018123006','gamma',6,2,[1:8],[]
            'TK18121401/2018123007','beta',6,7,[1:8],[]
            'TK18121401/2018123008','theta',6,7,[1:8],[]
            'TK18121401/2018123009','gamma',6,7,[1:8],[]
                 };
        
      case 'TK19011801',%%PD
        dataPath={
            'TK19011801/2019020601','beta',2,3,[1:8],[]
            'TK19011801/2019020602','theta',2,3,[1:8],[]
            'TK19011801/2019020603','gamma',2,3,[1:8],[]
            'TK19011801/2019020604','beta',2,8,[1:8],[]
            'TK19011801/2019020605','theta',2,8,[1:8],[]
            'TK19011801/2019020606','gamma',2,8,[1:8],[]
            'TK19011801/2019020704','beta',7,8,[1:8],[]
            'TK19011801/2019020705','theta',7,8,[1:8],[]
            'TK19011801/2019020706','gamma',7,8,[1:8],[]
            'TK19011801/2019020707','beta',7,3,[1:8],[]
            'TK19011801/2019020708','theta',7,3,[1:8],[]
            'TK19011801/2019020801','gamma',7,3,[1:8],[]
            'TK19011801/2019020802','theta',7,3,[1:8],[]
            'TK19011801/2019020803','beta',7,3,[1:8],[]
                 };

      case 'TK19012501',%%PD
        dataPath={
            'TK19012501/2019021204','beta',7,5,[1:8],[]
            'TK19012501/2019021205','theta',7,5,[1:8],[]
            'TK19012501/2019021206','gamma',7,5,[1:8],[]
            'TK19012501/2019021207','beta',7,1,[1:8],[]
            'TK19012501/2019021208','theta',7,1,[1:8],[]
            'TK19012501/2019021304','beta',2,1,[1:8],[]
            'TK19012501/2019021305','theta',2,1,[1:8],[]
            'TK19012501/2019021306','gamma',2,1,[1:8],[]
            'TK19012501/2019021401','beta',2,5,[1:8],[]
            'TK19012501/2019021402','theta',2,5,[1:8],[]
            'TK19012501/2019021403','gamma',2,5,[1:8],[]
            'TK19012501/2019021404','beta',2,1,[1:8],[]
            'TK19012501/2019021405','theta',2,1,[1:8],[]
            'TK19012501/2019021406','gamma',2,1,[1:8],[]
            'TK19012501/2019021407','gamma',7,1,[1:8],[]
            'TK19012501/2019021408','theta',7,1,[1:8],[]
            'TK19012501/2019021409','beta',7,1,[1:8],[]
                 };
    end

  case '2019',
    homePath='/mnt/TakahashiLabA/mouse/2019/';

    switch (an)
      case 'TK19040101',%%PD 
        dataPath={
            'TK19040101/2019041701','beta',6,1,[1:8],[]
            'TK19040101/2019041702','theta',6,1,[1:8],[]
            'TK19040101/2019041703','gamma',6,1,[1:8],[]
            'TK19040101/2019041704','beta',6,7,[1:8],[]
            'TK19040101/2019041705','theta',6,7,[1:8],[]
            'TK19040101/2019041706','gamma',6,7,[1:8],[]
            'TK19040101/2019041801','beta',3,2,[1:8],[]
            'TK19040101/2019041802','theta',3,2,[1:8],[]
            'TK19040101/2019041803','gamma',3,2,[1:8],[]
            'TK19040101/2019041804','beta',3,7,[1:8],[]
            %            'TK19040101/2019041805','theta',3,7,[1:8],[]
            'TK19040101/2019041806','gamma',3,7,[1:8],[]
            'TK19040101/2019041901','beta',3,3,[1:8],[]
            'TK19040101/2019041902','theta',3,3,[1:8],[]
            'TK19040101/2019041903','gamma',3,3,[1:8],[]
            'TK19040101/2019041904','gamma',3,7,[1:8],[]
            'TK19040101/2019041905','theta',3,7,[1:8],[]
            'TK19040101/2019041906','beta',3,7,[1:8],[]
            'TK19040101/2019042301','beta',6,1,[1:8],[]
            'TK19040101/2019042302','theta',6,1,[1:8],[]
            'TK19040101/2019042303','gamma',6,1,[1:8],[]
            'TK19040101/2019042304','beta',6,7,[1:8],[]
            'TK19040101/2019042305','theta',6,7,[1:8],[]
            'TK19040101/2019042306','gamma',6,7,[1:8],[]
                 };

      case 'TK19050801',%%PD 
        dataPath={
            'TK19050801/2019051404','beta',3,3,[1:8],[]
            'TK19050801/2019051405','theta',3,3,[1:8],[]
            'TK19050801/2019051406','gamma',3,3,[1:8],[]
            'TK19050801/2019051501','beta',6,4,[1:8],[]
            'TK19050801/2019051502','theta',6,4,[1:8],[]
            'TK19050801/2019051503','gamma',6,4,[1:8],[]
            'TK19050801/2019051504','gamma',6,7,[1:8],[]
            'TK19050801/2019051505','theta',6,7,[1:8],[]
            'TK19050801/2019051506','beta',6,7,[1:8],[]
            'TK19050801/2019051701','beta',3,4,[1:8],[]
            'TK19050801/2019051702','theta',3,4,[1:8],[]
            'TK19050801/2019051703','gamma',3,4,[1:8],[]
            'TK19050801/2019051704','gamma',3,7,[1:8],[]
            'TK19050801/2019051705','theta',3,7,[1:8],[]
            'TK19050801/2019051706','beta',3,7,[1:8],[]
            'TK19050801/2019052104','beta',6,4,[1:8],[]
            'TK19050801/2019052105','theta',6,4,[1:8],[]
            'TK19050801/2019052106','gamma',6,4,[1:8],[]
            'TK19050801/2019052107','gamma',6,5,[1:8],[]
            'TK19050801/2019052108','theta',6,5,[1:8],[]
            'TK19050801/2019052109','beta',6,5,[1:8],[]
                 };

      case 'TK19051601',%%PD
        dataPath={
            'TK19051601/2019052201','beta',3,3,[1:8],[]
            'TK19051601/2019052202','theta',3,3,[1:8],[]
            %            'TK19051601/2019052203','gamma',3,3,[1:8],[]
            'TK19051601/2019052204','beta',3,7,[1:8],[]
            'TK19051601/2019052205','theta',3,7,[1:8],[]
            %            'TK19051601/2019052206','gamma',3,7,[1:8],[]
            'TK19051601/2019052301','beta',6,3,[1:8],[]
            'TK19051601/2019052302','theta',6,3,[1:8],[]
            'TK19051601/2019052303','gamma',6,3,[1:8],[]
            'TK19051601/2019052304','beta',6,7,[1:8],[]
            'TK19051601/2019052305','theta',6,7,[1:8],[]
            'TK19051601/2019052306','gamma',6,7,[1:8],[]
            'TK19051601/2019052401','beta',3,2,[1:8],[]
            'TK19051601/2019052402','theta',3,2,[1:8],[]
            'TK19051601/2019052403','gamma',3,2,[1:8],[]
            'TK19051601/2019052404','beta',3,8,[1:8],[]
            'TK19051601/2019052405','theta',3,8,[1:8],[]
            'TK19051601/2019052406','gamma',3,8,[1:8],[]
                 };

      case 'TK19060301',%%PD 20220628
        dataPath={
            'TK19060301/2019061104','beta',6,6,[1:8],[]
            'TK19060301/2019061105','theta',6,6,[1:8],[]
            'TK19060301/2019061106','gamma',6,6,[1:8],[]
                 };

    end
    
    case '2019f',
        homePath='f:/2019/';

        switch (an)
          case 'TK19081901',%PD, no ChRs on Str
                dataPath={
                    'TK19081901/2019092701','beta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019092702','theta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019092703','gamma',6,8,[1:5 8],[6 7]
                    'TK19081901/2019092704','beta',6,1,[1:5 8],[6 7]
                    'TK19081901/2019092705','theta',6,1,[1:5 8],[6 7]
                    'TK19081901/2019092706','gamma',6,1,[1:5 8],[6 7]
                    'TK19081901/2019093001','beta',3,4,[1:5 8],[6 7]
                    'TK19081901/2019093002','theta',3,4,[1:5 8],[6 7]
                    'TK19081901/2019093003','gamma',3,4,[1:5 8],[6 7]
                    'TK19081901/2019093004','beta',3,8,[1:5 8],[6 7]
                    'TK19081901/2019093005','theta',3,8,[1:5 8],[6 7]
                    'TK19081901/2019093006','gamma',3,8,[1:5 8],[6 7]
                    'TK19081901/2019100101','beta',6,5,[1:5 8],[6 7]
                    'TK19081901/2019100102','theta',6,5,[1:5 8],[6 7]
                    'TK19081901/2019100103','gamma',6,5,[1:5 8],[6 7]
                    'TK19081901/2019100104','beta',6,4,[1:5 8],[6 7]
                    'TK19081901/2019100105','theta',6,4,[1:5 8],[6 7]
                    'TK19081901/2019100106','gamma',6,4,[1:5 8],[6 7]
                    'TK19081901/2019100301','beta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100302','theta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100303','gamma',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100304','beta',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100305','theta',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100306','gamma',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100701','beta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100702','theta',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100703','gamma',6,8,[1:5 8],[6 7]
                    'TK19081901/2019100704','beta',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100705','theta',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100706','gamma',6,2,[1:5 8],[6 7]
                    'TK19081901/2019100801','beta',3,4,[1:5 7 8],[6]
                    'TK19081901/2019100802','theta',3,4,[1:5 7 8],[6]
                    'TK19081901/2019100803','gamma',3,4,[1:5 7 8],[6]
                    'TK19081901/2019100804','beta',3,8,[1:5 7 8],[6]
                    'TK19081901/2019100805','theta',3,8,[1:5 7 8],[6]
                    'TK19081901/2019100806','gamma',3,8,[1:5 7 8],[6]
                         };

          case 'TK19100801',%PD,Str-Str(1 8), Ctx-LED
                dataPath={
                    'TK19100801/2019102401','beta',7,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102402','theta',7,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102403','gamma',7,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102501','beta',7,4,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102502','theta',7,4,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102503','gamma',7,4,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102801','beta',2,1,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102802','theta',2,1,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102803','gamma',2,1,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102804','beta',2,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102805','theta',2,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019102806','gamma',2,6,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019103001','beta',7,8,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019103002','theta',7,8,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019103003','gamma',7,8,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019103101','beta',7,3,[3 4 5 6 7],[1 2 8]
                    %                    'TK19100801/2019103102','theta',7,3,[3 4 5 6 7],[1 2 8]
                    %'TK19100801/2019103103','gamma',7,3,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110101','beta',2,3,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110102','theta',2,3,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110103','gamma',2,3,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110501','beta',2,8,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110502','theta',2,8,[3 4 5 6 7],[1 2 8]
                    'TK19100801/2019110503','gamma',2,8,[3 4 5 6 7],[1 2 8]
                         };

        case 'AZ19091201',%normal, ChR2
                dataPath={
                    'AZ19091201/2019101001','beta',6,7,[1:8],[]
                    'AZ19091201/2019101002','theta',6,7,[1:8],[]
                    'AZ19091201/2019101003','gamma',6,7,[1:8],[]
                    'AZ19091201/2019101101','beta',6,1,[1:8],[]
                    'AZ19091201/2019101102','theta',6,1,[1:8],[]
                    'AZ19091201/2019101103','gamma',6,1,[1:8],[]
                    'AZ19091201/2019101501','beta',3,1,[1:8],[]
                    'AZ19091201/2019101502','theta',3,1,[1:8],[]
                    'AZ19091201/2019101503','gamma',3,1,[1:8],[]
                    'AZ19091201/2019101601','beta',3,5,[1:8],[]
                    'AZ19091201/2019101602','theta',3,5,[1:8],[]
                    'AZ19091201/2019101603','gamma',3,5,[1:8],[]
                         };

          case 'AZ19120501',%PD
                dataPath={
                    'AZ19120501/2020012101','beta',3,1,[1:8],[]
                    'AZ19120501/2020012102','theta',3,1,[1:8],[]
                    'AZ19120501/2020012103','gamma',3,1,[1:8],[]
                    'AZ19120501/2020012201','beta',3,7,[1:8],[]
                    'AZ19120501/2020012202','theta',3,7,[1:8],[]
                    'AZ19120501/2020012203','gamma',3,7,[1:8],[]
                    'AZ19120501/2020012301','beta',6,7,[1:8],[]
                    'AZ19120501/2020012302','theta',6,7,[1:8],[]
                    'AZ19120501/2020012303','gamma',6,7,[1:8],[]
                    'AZ19120501/2020012401','beta',6,1,[1:8],[]
                    'AZ19120501/2020012402','theta',6,1,[1:8],[]
                    'AZ19120501/2020012403','gamma',6,1,[1:8],[]
                    'AZ19120501/2020012801','beta',3,1,[1:8],[]
                    'AZ19120501/2020012802','theta',3,1,[1:8],[]
                    'AZ19120501/2020012803','gamma',3,1,[1:8],[]
                    'AZ19120501/2020020301','beta',3,7,[1:8],[]
                    'AZ19120501/2020020302','theta',3,7,[1:8],[]
                    'AZ19120501/2020020303','gamma',3,7,[1:8],[]
                    'AZ19120501/2020020401','beta',6,7,[1:8],[]
                    'AZ19120501/2020020402','theta',6,7,[1:8],[]
                    'AZ19120501/2020020403','gamma',6,7,[1:8],[]
                    'AZ19120501/2020020501','beta',6,1,[1:8],[]
                    'AZ19120501/2020020502','theta',6,1,[1:8],[]
                    'AZ19120501/2020020503','gamma',6,1,[1:8],[]
                         };

          case 'AZ20013101',%PD
                dataPath={
                    'AZ20013101/2020021702','beta',6,8,[1:8],[]
                    'AZ20013101/2020021703','theta',6,8,[1:8],[]
                    'AZ20013101/2020021704','gamma',6,8,[1:8],[]
                    'AZ20013101/2020021801','beta',6,1,[1:8],[]
                    'AZ20013101/2020021802','theta',6,1,[1:8],[]
                    'AZ20013101/2020021803','gamma',6,1,[1:8],[]
                         };

          case 'AZ20031201',%PD
                dataPath={
                    'AZ20031201/2020032301','beta',3,4,[1:8],[]
                    'AZ20031201/2020032302','theta',3,4,[1:8],[]
                    'AZ20031201/2020032303','gamma',3,4,[1:8],[]
                    'AZ20031201/2020032401','beta',3,2,[1:8],[]
                    'AZ20031201/2020032402','theta',3,2,[1:8],[]
                    'AZ20031201/2020032403','gamma',3,2,[1:8],[]
                    'AZ20031201/2020032404','beta',3,5,[1:8],[]
                    'AZ20031201/2020032405','theta',3,5,[1:8],[]
                    'AZ20031201/2020032406','gamma',3,5,[1:8],[]
                    'AZ20031201/2020032501','beta',6,8,[1:8],[]
                    'AZ20031201/2020032502','theta',6,8,[1:8],[]
                    'AZ20031201/2020032503','gamma',6,8,[1:8],[]
                    'AZ20031201/2020032504','beta',6,3,[1:8],[]
                    'AZ20031201/2020032601','beta',6,1,[1:8],[]
                    'AZ20031201/2020032602','theta',6,1,[1:8],[]
                    'AZ20031201/2020032603','gamma',6,1,[1:8],[]
                    'AZ20031201/2020033001','beta',3,3,[1:8],[]
                    'AZ20031201/2020033002','theta',3,3,[1:8],[]
                    'AZ20031201/2020033003','gamma',3,3,[1:8],[]
                    'AZ20031201/2020033004','beta',3,7,[1:8],[]
                    'AZ20031201/2020033101','beta',3,6,[1:8],[]
                    'AZ20031201/2020033102','theta',3,6,[1:8],[]
                    'AZ20031201/2020033103','gamma',3,6,[1:8],[]
                    'AZ20031201/2020040101','beta',6,7,[1:8],[]
                    'AZ20031201/2020040102','theta',6,7,[1:8],[]
                    'AZ20031201/2020040103','gamma',6,7,[1:8],[]
                    'AZ20031201/2020040104','beta',6,4,[1:8],[]
                    'AZ20031201/2020040105','theta',6,4,[1:8],[]
                    'AZ20031201/2020040106','gamma',6,4,[1:8],[]
                    %
                    };

        end
    
  case '2020h',
        homePath='h:/2020/';

        switch (an)        
          case 'AZ20030601',%PD
                dataPath={
                    'AZ20030601/2020040601','beta',3,2,[1:8],[]
                    'AZ20030601/2020040602','theta',3,2,[1:8],[]
                    'AZ20030601/2020040603','gamma',3,2,[1:8],[]
                    'AZ20030601/2020040604','beta',3,8,[1:8],[]
                    'AZ20030601/2020040605','theta',3,8,[1:8],[]
                    'AZ20030601/2020040606','gamma',3,8,[1:8],[]
                    'AZ20030601/2020040701','beta',6,7,[1:8],[]
                    'AZ20030601/2020040702','theta',6,7,[1:8],[]
                    'AZ20030601/2020040703','gamma',6,7,[1:8],[]
                    'AZ20030601/2020040704','beta',6,3,[1:8],[]
                    'AZ20030601/2020040705','theta',6,3,[1:8],[]
                    'AZ20030601/2020040706','gamma',6,3,[1:8],[]
                    };

          case 'AZ20031901',%PD
                dataPath={
                    'AZ20031901/2020040901','beta',3,2,[1:8],[]
                    'AZ20031901/2020040902','theta',3,2,[1:8],[]
                    'AZ20031901/2020040903','gamma',3,2,[1:8],[]
                    'AZ20031901/2020040904','beta',3,7,[1:8],[]
                    'AZ20031901/2020040905','theta',3,7,[1:8],[]
                    'AZ20031901/2020040906','gamma',3,7,[1:8],[]
                    'AZ20031901/2020041001','beta',6,8,[1:8],[]
                    'AZ20031901/2020041002','theta',6,8,[1:8],[]
                    'AZ20031901/2020041003','gamma',6,8,[1:8],[]
                    'AZ20031901/2020041004','beta',6,3,[1:8],[]
                    'AZ20031901/2020041005','theta',6,3,[1:8],[]
                    'AZ20031901/2020041006','gamma',6,3,[1:8],[]
                         };
                
          case 'AZ20070101',%PD
                dataPath={
                    'AZ20070101/2020072201','beta',3,2,[1:3 7 8],[4:6]
                    'AZ20070101/2020072202','theta',3,2,[1:3 7 8],[4:6]
                    'AZ20070101/2020072203','gamma',3,2,[1:3 7 8],[4:6]
                    'AZ20070101/2020072204','beta',3,8,[1:3 7 8],[4:6]
                    'AZ20070101/2020072205','theta',3,8,[1:3 7 8],[4:6]
                    'AZ20070101/2020072206','gamma',3,8,[1:3 7 8],[4:6]
                    'AZ20070101/2020072801','beta',6,8,[1:3 7 8],[4:6] %nostim
                    'AZ20070101/2020072802','theta',6,8,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020072803','gamma',6,8,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020072804','beta',6,3,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020072805','theta',6,3,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020072806','gamma',6,3,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020073001','beta',3,4,[1:3 7 8],[4:6]
                    'AZ20070101/2020073002','theta',3,4,[1:3 7 8],[4:6]
                    'AZ20070101/2020073003','gamma',3,4,[1:3 7 8],[4:6]
                    'AZ20070101/2020073004','beta',3,6,[1:3 7 8],[4:6]
                    'AZ20070101/2020073005','theta',3,6,[1:3 7 8],[4:6]
                    'AZ20070101/2020073006','gamma',3,6,[1:3 7 8],[4:6]
                    'AZ20070101/2020073101','beta',6,6,[1:3 7 8],[4:6]%nostim
                    'AZ20070101/2020073102','theta',6,6,[1:3 7 8],[4:6]%nostim
                    %                    'AZ20070101/2020073103','gamma',6,6,[1:3 7 8],[4:6]
                    'AZ20070101/2020073104','beta',6,2,[1:3 7 8],[4:6]
                    'AZ20070101/2020073105','theta',6,2,[1:3 7 8],[4:6]
                    'AZ20070101/2020073106','gamma',6,2,[1:3 7 8],[4:6]
                    };
        end

  case '2020e',
        homePath='e:/2020/';

        switch (an)        
          case 'AZ20121501',%PD
                dataPath={
                    'AZ20121501/2021011501','beta',3,3,[1:8],[]
                    'AZ20121501/2021011502','theta',3,3,[1:8],[]
                    'AZ20121501/2021011503','gamma',3,3,[1:8],[]
                    'AZ20121501/2021011504','beta',3,5,[1:8],[]
                    'AZ20121501/2021011505','theta',3,5,[1:8],[]
                    'AZ20121501/2021011506','gamma',3,5,[1:8],[]
                    'AZ20121501/2021011801','beta',6,5,[1:8],[]
                    'AZ20121501/2021011802','theta',6,5,[1:8],[]
                    'AZ20121501/2021011803','gamma',6,5,[1:8],[]
                    'AZ20121501/2021011804','beta',6,3,[1:8],[]
                    'AZ20121501/2021011805','theta',6,3,[1:8],[]
                    'AZ20121501/2021011806','gamma',6,3,[1:8],[]
                    'AZ20121501/2021012101','beta',3,1,[1:8],[]
                    'AZ20121501/2021012102','theta',3,1,[1:8],[]
                    'AZ20121501/2021012103','gamma',3,1,[1:8],[]
                    'AZ20121501/2021012104','beta',3,6,[1:8],[]
                    'AZ20121501/2021012105','theta',3,6,[1:8],[]
                    'AZ20121501/2021012106','gamma',3,6,[1:8],[]
                    'AZ20121501/2021012201','beta',6,7,[1:8],[]
                    'AZ20121501/2021012202','theta',6,7,[1:8],[]
                    'AZ20121501/2021012601','beta',6,1,[1:8],[]
                    'AZ20121501/2021012602','theta',6,1,[1:8],[]
                    'AZ20121501/2021012603','gamma',6,1,[1:8],[]
                    };

          case 'AZ21012601',%PD
                dataPath={
                    'AZ21012601/2021020401','beta',3,2,[1:5 7:8],[6]
                    'AZ21012601/2021020402','theta',3,2,[1:5 7:8],[6]
                    'AZ21012601/2021020403','gamma',3,2,[1:5 7:8],[6]
                    'AZ21012601/2021020404','beta',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021020405','theta',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021020406','gamma',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021020501','beta',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021020502','theta',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021020503','gamma',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021020504','beta',6,3,[1:5 7:8],[6]
                    'AZ21012601/2021020505','theta',6,3,[1:5 7:8],[6]
                    'AZ21012601/2021020506','gamma',6,3,[1:5 7:8],[6]
                    'AZ21012601/2021020801','beta',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021020802','theta',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021020803','gamma',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021020804','beta',3,6,[1:5 7:8],[6]
                    'AZ21012601/2021020805','theta',3,6,[1:5 7:8],[6]
                    'AZ21012601/2021020806','gamma',3,6,[1:5 7:8],[6]
                    'AZ21012601/2021020901','beta',6,5,[1:5 7:8],[6]
                    'AZ21012601/2021020902','theta',6,5,[1:5 7:8],[6]
                    'AZ21012601/2021020903','gamma',6,5,[1:5 7:8],[6]
                    'AZ21012601/2021020904','beta',6,1,[1:5 7:8],[6]
                    'AZ21012601/2021020905','theta',6,1,[1:5 7:8],[6]
                    'AZ21012601/2021020906','gamma',6,1,[1:5 7:8],[6]
                    'AZ21012601/2021021501','beta',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021021502','theta',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021021503','gamma',3,4,[1:5 7:8],[6]
                    'AZ21012601/2021021504','beta',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021021505','theta',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021021506','gamma',3,7,[1:5 7:8],[6]
                    'AZ21012601/2021021601','beta',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021021602','theta',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021021603','gamma',6,8,[1:5 7:8],[6]
                    'AZ21012601/2021021604','beta',6,1,[1:5 7:8],[6]
                    'AZ21012601/2021021605','theta',6,1,[1:5 7:8],[6]
                    'AZ21012601/2021021606','gamma',6,1,[1:5 7:8],[6]
                         };

          case 'AZ21010701',%PD
                dataPath={
                    'AZ21010701/2021030201','beta',3,2,[1:8],[]
                    'AZ21010701/2021030202','theta',3,2,[1:8],[]
                    'AZ21010701/2021030203','gamma',3,2,[1:8],[]
                    'AZ21010701/2021030204','beta',3,5,[1:8],[]
                    'AZ21010701/2021030205','theta',3,5,[1:8],[]
                    'AZ21010701/2021030206','gamma',3,5,[1:8],[]
                    'AZ21010701/2021030301','beta',6,5,[1:8],[]
                    'AZ21010701/2021030302','theta',6,5,[1:8],[]
                    'AZ21010701/2021030303','gamma',6,5,[1:8],[]
                    'AZ21010701/2021030401','beta',6,5,[1:8],[]
                    'AZ21010701/2021030402','theta',6,5,[1:8],[]
                    'AZ21010701/2021030403','gamma',6,5,[1:8],[]
                    'AZ21010701/2021030404','beta',6,2,[1:8],[]
                    'AZ21010701/2021030405','theta',6,2,[1:8],[]
                    'AZ21010701/2021030406','gamma',6,2,[1:8],[]
                         };

          case 'AZ21012501',%PD
                dataPath={
                %                    'AZ21012501/2021031001','beta',3,1,[1:2 4 5 7],[3 6 8]
                    'AZ21012501/2021031002','theta',3,1,[1:8],[]
                    'AZ21012501/2021031003','gamma',3,1,[1:8],[]
                    'AZ21012501/2021031004','beta',3,7,[1:8],[]
                    'AZ21012501/2021031005','theta',3,7,[1:8],[]
                    'AZ21012501/2021031006','gamma',3,7,[1:8],[]
                         };
        end


  case '2022h',
        homePath='h:/2022/';

        switch (an)        
          case 'AZ22070601',%normal control
                dataPath={
                    'AZ22070601/2022072601','beta',3,3,[1:8],[]
                    'AZ22070601/2022072602','theta',3,3,[1:8],[]
                    'AZ22070601/2022072603','gamma',3,3,[1:8],[]
                    'AZ22070601/2022072604','beta',3,5,[1:8],[]
                    'AZ22070601/2022072606','theta',3,5,[1:8],[]
                    'AZ22070601/2022072607','gamma',3,5,[1:8],[]
                    'AZ22070601/2022072701','beta',6,7,[1:8],[]
                    'AZ22070601/2022072702','theta',6,7,[1:8],[]
                    'AZ22070601/2022072703','gamma',6,7,[1:8],[]
                    'AZ22070601/2022072704','beta',6,4,[1:8],[]
                    'AZ22070601/2022072705','theta',6,4,[1:8],[]
                    'AZ22070601/2022072706','gamma',6,4,[1:8],[]
                         };
          case 'AZ22070701',%normal control
                dataPath={
                    'AZ22070701/2022072901','beta',3,1,[1:8],[]
                    'AZ22070701/2022072902','theta',3,1,[1:8],[]
                    'AZ22070701/2022072903','gamma',3,1,[1:8],[]
                    'AZ22070701/2022072904','beta',3,8,[1:8],[]
                    'AZ22070701/2022072905','theta',3,8,[1:8],[]
                    'AZ22070701/2022072906','gamma',3,8,[1:8],[]
                    'AZ22070701/2022080101','beta',6,8,[1:8],[]
                    'AZ22070701/2022080102','theta',6,8,[1:8],[]
                    'AZ22070701/2022080103','gamma',6,8,[1:8],[]
                    'AZ22070701/2022080105','beta',6,1,[1:8],[]
                    'AZ22070701/2022080106','theta',6,1,[1:8],[]
                    'AZ22070701/2022080107','gamma',6,1,[1:8],[]
                         };
        end

  case '2022g',
        homePath='g:/2022/';

        switch (an)        
          case 'AZ22080501',%normal control Ctx
                dataPath={
                    'AZ22080501/2022082501','beta',3,1,[1:8],[]
                    'AZ22080501/2022082502','theta',3,1,[1:8],[]
                    'AZ22080501/2022082503','gamma',3,1,[1:8],[]
                    'AZ22080501/2022082504','beta',3,6,[1:8],[]
                    'AZ22080501/2022082505','theta',3,6,[1:8],[]
                    'AZ22080501/2022082506','gamma',3,6,[1:8],[]
                    'AZ22080501/2022082901','beta',6,7,[1:8],[]
                    'AZ22080501/2022082902','theta',6,7,[1:8],[]
                    'AZ22080501/2022082903','gamma',6,7,[1:8],[]
                    'AZ22080501/2022082904','beta',6,1,[1:8],[]
                    'AZ22080501/2022082905','theta',6,1,[1:8],[]
                    'AZ22080501/2022082906','gamma',6,1,[1:8],[]
                         };
case 'AZ22091401',%normal control,ChR2
                dataPath={
                    'AZ22091401/2022100601','beta',3,2,[1:8],[]
                    'AZ22091401/2022100602','theta',3,2,[1:8],[]
                    'AZ22091401/2022100603','gamma',3,2,[1:8],[]
                    'AZ22091401/2022100604','beta',3,7,[1:8],[]
                    'AZ22091401/2022100605','theta',3,7,[1:8],[]
                    'AZ22091401/2022100606','gamma',3,7,[1:8],[]
                    'AZ22091401/2022100701','beta',6,7,[1:8],[]
                    'AZ22091401/2022100702','theta',6,7,[1:8],[]
                    'AZ22091401/2022100703','gamma',6,7,[1:8],[]
                    'AZ22091401/2022100704','beta',6,1,[1:8],[]
                    'AZ22091401/2022100705','theta',6,1,[1:8],[]
                    'AZ22091401/2022100706','gamma',6,1,[1:8],[]
                         };


case 'AZ22091402',%normal control,ChR2
                dataPath={
                    'AZ22091402/2022101801','beta',3,3,[1:8],[]
                    'AZ22091402/2022101802','theta',3,3,[1:8],[]
                    'AZ22091402/2022101803','gamma',3,3,[1:8],[]
                    'AZ22091402/2022101804','beta',3,5,[1:8],[]
                    'AZ22091402/2022101805','theta',3,5,[1:8],[]
                    'AZ22091402/2022101806','gamma',3,5,[1:8],[]
                    'AZ22091402/2022101901','beta',6,5,[1:8],[]
                    'AZ22091402/2022101902','theta',6,5,[1:8],[]
                    'AZ22091402/2022101903','gamma',6,5,[1:8],[]
                    'AZ22091402/2022101904','beta',6,3,[1:8],[]
                    'AZ22091402/2022101905','theta',6,3,[1:8],[]
                    'AZ22091402/2022101906','gamma',6,3,[1:8],[]
                         };
        end

end


