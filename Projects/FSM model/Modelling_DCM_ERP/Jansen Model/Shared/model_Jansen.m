function dxdt=model_Jansen(t,x)

global Ae2 Ai2 a2 DD a1 b2 b1 Ae Ai ae ai c1 c2 c3 c4 c5 c6 c7 c8 c9 th P S A c sel Neq


dxdt = zeros(Neq,1);
DD = 0;
if sel==1
    
    %% without feedback
    dxdt(1) = x(4);
    dxdt(2) = x(5);
    dxdt(3) = x(6);
    dxdt(4) = Ae*ae*(P + c2*S(c1*x(3)))- 2*ae*x(4) - ae^2*x(1);
    dxdt(5) = Ai*ai*c4*S(c3*x(3))      - 2*ai*x(5) - ai^2*x(2);
    dxdt(6) = Ae*ae*S(x(1)-x(2))       - 2*ae*x(6) - ae^2*x(3);
    dxdt(7) = x(4)-x(5);
    
elseif sel==2
    
    
    %% Feedback in In-in
    dxdt(1) = x(6);
    dxdt(2) = x(7);
    dxdt(3) = x(8);
    dxdt(4) = x(9);
    dxdt(5) = x(10);
    
    dxdt(6) = Ae*ae*(c1*S(x(2)-x(3)-DD)+P)- 2*ae*x(6) - ae^2*x(1);
    dxdt(7) = Ae*ae*c2*S(x(1))         - 2*ae*x(7) - ae^2*x(2);
    dxdt(8) = Ai*ai*c3*S(x(4)-x(5))    - 2*ai*x(8) - ai^2*x(3);
    dxdt(9) = Ae*ae*c4*S(x(2)-x(3))    - 2*ae*x(9) - ae^2*x(4);
    dxdt(10)= Ai*ai*c5*S(x(4)-x(5))    - 2*ai*x(10)- ai^2*x(5);
    dxdt(11)= x(7)-x(8);
    
elseif sel==3
    
    %% Feedback in Excitary pyrimidal
    %     y1 = x(3)-x(5);
    %     y2 = x(3);
    %     y3 = x(1)-x(2);
    %
    dxdt(1) = x(6);
    dxdt(2) = x(7);
    dxdt(3) = x(8);
    dxdt(4) = x(9);
    dxdt(5) = x(10);
    
    dxdt(6) = Ae*ae*(P+c2*S(c1*y1))-2*ae*x(6)-ae^2*x(1);
    dxdt(7) = Ai*ai*(c4*S(c3*y2))-2*ai*x(7)-ai^2*x(2);
    dxdt(8) = Ae*ae*S(y3)-2*ae*x(8)-ae^2*x(3);
    dxdt(9) = Ai*ai*(c5*S(c3*y2))-2*ai*x(9) -ai^2*x(4);
    dxdt(10)= Ae*ae*(c6*S(c1*y1))-2*ae*x(10)-ae^2*x(5);
    
elseif sel==4
    
    %% Feedback loop in In-in and Ex-py
    
    dxdt(1) = x(7);
    dxdt(2) = x(8);
    dxdt(3) = x(9);
    dxdt(4) = x(10);
    dxdt(5) = x(11);
    dxdt(6) = x(12);
    
    dxdt(7) = Ae*ae*(c1*S(x(2)-x(3))+P)- 2*ae*x(7)  - ae^2*x(1);
    dxdt(8) = Ae*ae* c2*S(x(1)-x(6))   - 2*ae*x(8)  - ae^2*x(2);
    dxdt(9) = Ai*ai* c3*S(x(4)-x(5))   - 2*ai*x(9)  - ai^2*x(3);
    dxdt(10)= Ae*ae* c4*S(x(2)-x(3))   - 2*ae*x(10) - ae^2*x(4);
    dxdt(11)= Ai*ai* c5*S(x(4)-x(5))   - 2*ai*x(11) - ai^2*x(5);
    dxdt(12)= Ae*ae* c6*S(x(1)-x(6))   - 2*ae*x(12) - ae^2*x(6);
    dxdt(13)= x(8)-x(9);
    
    
elseif sel==5
    
    %% Feedback loop in In-in and Ex-py and E-in
    dxdt(1) = x(8);
    dxdt(2) = x(9);
    dxdt(3) = x(10);
    dxdt(4) = x(11);
    dxdt(5) = x(12);
    dxdt(6) = x(13);
    dxdt(7) = x(14);
    
    dxdt(8) = Ae*ae*(c1*S(x(2)-x(3)+x(7)-DD)+P)  - 2*ae*x(8)  - ae^2*x(1);
    dxdt(9) = Ae*ae* c2*S(x(1)+x(6))     - 2*ae*x(9)  - ae^2*x(2);
    dxdt(10)= Ai*ai* c3*S(x(4)-x(5))     - 2*ai*x(10) - ai^2*x(3);
    dxdt(11)= Ae*ae* c4*S(x(2)-x(3)+x(7)-DD)- 2*ae*x(11) - ae^2*x(4);
    dxdt(12)= Ai*ai* c5*S(x(4)-x(5))     - 2*ai*x(12) - ai^2*x(5);
    dxdt(13)= Ae*ae* c6*S(x(1)+x(6))     - 2*ae*x(13) - ae^2*x(6);
    dxdt(14)= Ae*ae* c7*S(x(2)-x(3)+x(7)-DD)- 2*ae*x(14) - ae^2*x(7);
    dxdt(15)= x(9)-x(10)+x(14);
    
elseif sel==6
    
    %% Feedback loop in In-in and Ex-py and E-in
    dxdt(1) = x(8);
    dxdt(2) = x(9);
    dxdt(3) = x(10);
    dxdt(4) = x(11);
    dxdt(5) = x(12);
    dxdt(6) = x(13);
    dxdt(7) = x(14);
    
    dxdt(8) = Ae2*(a2-a1)*(c1*S(x(2)-x(3)+x(7))+P)  - 2*(a2-a2)*x(8)  - (a1*a2)*x(1);
    dxdt(9) = Ae2*(a2-a1)* c2*S(x(1)+x(6))          - 2*(a2-a1)*x(9)  - (a1*a2)*x(2);
    dxdt(10)= Ai2*(b2-b1)* c3*S(x(4)-x(5))          - 2*(b2-b1)*x(10) - (b1*b2)*x(3);
    dxdt(11)= Ae2*(a2-a1)* c4*S(x(2)-x(3)+x(7))     - 2*(a2-a1)*x(11) - (a1*a2)*x(4);
    dxdt(12)= Ai2*(b2-b1)* c5*S(x(4)-x(5))          - 2*(b2-b1)*x(12) - (b1*b2)*x(5);
    dxdt(13)= Ae2*(a2-a1)* c6*S(x(1)+x(6))          - 2*(a2-a1)*x(13) - (a1*a2)*x(6);
    dxdt(14)= Ae2*(a2-a1)* c7*S(x(2)-x(3)+x(7))     - 2*(a2-a1)*x(14) - (a1*a2)*x(7);
    dxdt(15)= x(9)-x(10)+x(14);
    
    
elseif sel==7
    
    %% Feedback loop in In-in and Ex-py and E-in
    dxdt(1) = x(8);
    dxdt(8) = (He.*((A{1} + A{3})*S(x(15)) + c(1).*(S(x(2)-x(3)))+U)    - 2*x(8)  - x(1)./Te)./Te;
    %
    dxdt(2) = x(9);
    dxdt(9) = (He.*((A{2} + A{3})*S(x(15)) + c(2).*(S(x(1)-x(6))))      - 2*x(9)  - x(2)./Te)./Te;
    %
    dxdt(3) = x(10);
    dxdt(10)= (Hi*c(3).*(S(x(4)-x(5)))    - 2*x(10) - x(3)./Ti)./Ti;
    %
    dxdt(4) = x(11);
    dxdt(11)= (He.*((A{2} + A{3})*S(x(15))+ c(4)*S(x(2)-x(3)-x(7))- 2*x(11) - x(4)))./Te;
    %
    dxdt(5) = x(12);
    dxdt(12)= (Hi*c(5).*(S(x(4)-x(5)))         - 2*x(12) - x(5)./Ti)./Ti;
    %
    dxdt(6) = x(13);
    dxdt(13)= (He.*c(6).*(S(x(1)-x(6)))        - 2*x(13) - x(6)./Te)./Te;
    %
    dxdt(7) = x(14);
    dxdt(14)= (He.*c(7).*(S(x(2)-x(3)-x(7)))   - 2*x(14) - x(7)./Te)./Te;
    %
    dxdt(15)= x(9)-x(10)+x(14);
    
elseif sel==9
    
    
    dxdt(1) = x(8);
    dxdt(2) = x(9);
    dxdt(3) = x(10);
    dxdt(4) = x(11);
    dxdt(5) = x(12);
    dxdt(6) = x(13);
    dxdt(7) = x(14);
    
    dxdt(8) = Ae*ae*(c1*S(x(2)-x(3)+x(7)+x(22))+ P)  - 2*ae*x(8)  - ae^2*x(1);
    dxdt(9) = Ae*ae* c2*S(x(1)+x(6))     - 2*ae*x(9)  - ae^2*x(2);
    dxdt(10)= Ai*ai* c3*S(x(4)-x(5))     - 2*ai*x(10) - ai^2*x(3);
    dxdt(11)= Ae*ae* c4*S(x(2)-x(3)+x(7))- 2*ae*x(11) - ae^2*x(4);
    dxdt(12)= Ai*ai* c5*S(x(4)-x(5))     - 2*ai*x(12) - ai^2*x(5);
    dxdt(13)= Ae*ae* c6*S(x(1)+x(6))     - 2*ae*x(13) - ae^2*x(6);
    dxdt(14)= Ae*ae* c7*S(x(2)-x(3)+x(7))- 2*ae*x(14) - ae^2*x(7);
    dxdt(15)= x(9)-x(10)+x(14);
    
    dxdt(16)= x(17);
    dxdt(17)= Ae*ae* th*S(x(15))- 2*ae*x(17) - ae^2*x(16);
    dxdt(18)= x(19);
    dxdt(19)= Ai*ai*c9*S(c8*x(21))- 2*ai*x(19) - ai^2*x(18);
    dxdt(20)= x(21);
    dxdt(21)= Ae*ae* S(x(16)-x(18))- 2*ae*x(21) - ae^2*x(20);
    dxdt(22)= x(16)-x(18);
    
end

end