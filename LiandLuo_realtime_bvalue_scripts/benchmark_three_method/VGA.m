% The visibility degree derived by the “divide & conquer” algorithm. The strategy 
% utilizes the fact that the nodes on one side of the maximum cannot “see” another
% node on the other side. The algorithm finds the connections of the maximum node 
% with other nodes, and then divide the time series into two separated parts based
% on the maximum node, and repeat the above procedure until every segment contains
% only two nodes.

function [VG]=VGA(Dt,Dm,left,right,VG)
    [max_value,max_place0]=max(Dm(left:right));
    max_place=max_place0+left-1;
     for i=left:max_place-2
         tem_Dm=Dm(i+1:max_place-1);
         tem_Dt=Dt(i+1:max_place-1);
         cri_Dm=Dm(i)+(Dm(max_place)-Dm(i))*(tem_Dt-Dt(i))/(Dt(max_place)-Dt(i));
         cri=tem_Dm-cri_Dm;
         if max(cri)<0
            VG(i)=VG(i)+1;
            VG(max_place)=VG(max_place)+1;
         end
     end
     for i=max_place+2:right
         tem_Dm=Dm(max_place+1:i-1);
         tem_Dt=Dt(max_place+1:i-1);
         cri_Dm=Dm(i)+(Dm(max_place)-Dm(i))*(tem_Dt-Dt(i))/(Dt(max_place)-Dt(i));
         cri=tem_Dm-cri_Dm;
         if max(cri)<0
            VG(i)=VG(i)+1;
            VG(max_place)=VG(max_place)+1;
         end
     end
     if left<max_place-2
         VG=VGA(Dt,Dm,left,max_place-1,VG);
     end
     if right>max_place+2
         VG=VGA(Dt,Dm,max_place+1,right,VG);
     end
end