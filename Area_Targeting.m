function Area_Targeting(A1,A2)

% TotalArea only works for BCC, it calculates the approximate heat
% exchanger area required

% input A1 for hot streams in the format
%   [CPstream1 Ts Tt HTC; CPstream2 Ts Tt HTC;...]

% input A2 for cold streams in the format
%   [CPstream1 Ts Tt HTC; CPstream2 Ts Tt HTC;...]

temperatureinterval = forArea(A1,A2);

A = 0; % Area
a = 0; % lower T for Hot CC
t = 0; % lower T for Cold CC

for i = 1:numel(temperatureinterval)
    [~, xhot, yhot] = HotCC(A1); % get x and y coordinate values for hot CC
    [~, xcold,ycold] = ColdCC(A2); % get x and y coordinate values for cold CC
    
    TlowerforHotCC = a;
    TlowerforColdCC = t;
    
    % get Tupper for hot and cold CC
    [~,TupperforHotCC] = polyxpoly(xhot,yhot,[temperatureinterval(i)...
                                            temperatureinterval(i)],[0 1000]);
                                        
    [~,TupperforColdCC] = polyxpoly(xcold,ycold,[temperatureinterval(i)...
                                            temperatureinterval(i)],[0 1000]);
                                        
    a = max(TupperforHotCC);
    t = max(TupperforColdCC);
    
    
    TupperforHotCC = min(TupperforHotCC);
    TupperforColdCC = min(TupperforColdCC);
    
    if i == 1 % calculating area for the 1st segment
        [~,TlowerforHotCC] = polyxpoly(xhot,yhot,[0 0],[0 1000]);
        [~,TlowerforColdCC] = polyxpoly(xcold,ycold,[0 0],[0 1000]);
        
        % calculating log mean T
        LMT = ((TupperforHotCC - TupperforColdCC)-(TlowerforHotCC - TlowerforColdCC))...
            /log((TupperforHotCC - TupperforColdCC)/(TlowerforHotCC - TlowerforColdCC));
        
        % calculating Q/U in the formula A = Q/U * 1/LMT
        sumforhot = 0;
        sumforcold = 0;
        for j = 1:numel(A1(:,1))
            if A1(j,2) >= TupperforHotCC && A1(j,3) <= TlowerforHotCC
                sumforhot=sumforhot+A1(j,1)*(TupperforHotCC-TlowerforHotCC)/A1(j,4);
            end
        end
        
        for k = 1:numel(A2(:,1))
            if A2(k,2) <= TlowerforColdCC && A2(k,3) >= TupperforColdCC
                sumforcold=sumforcold+A2(k,1)*(TupperforColdCC-TlowerforColdCC)/A2(k,4);
            end
        end
        
        % calculating area using the formula A = Q/U * 1/LMT
        A = A + 1/LMT * (sumforhot + sumforcold);
        
        
    elseif i == numel(temperatureinterval) % calculating area for last segment
        
        % calculating log mean T
        LMT = ((TupperforHotCC - TupperforColdCC)-(TlowerforHotCC - TlowerforColdCC))...
            /log((TupperforHotCC - TupperforColdCC)/(TlowerforHotCC - TlowerforColdCC));
        
        % calculating Q/U in the formula A = Q/U * 1/LMT
        sumforhot = 0;
        sumforcold = 0;
        for j = 1:numel(A1(:,1))
            if A1(j,2) >= TupperforHotCC && A1(j,3) <= TlowerforHotCC
                sumforhot=sumforhot+A1(j,1)*(TupperforHotCC-TlowerforHotCC)/A1(j,4);
            end
        end
        
        for k = 1:numel(A2(:,1))
            if A2(k,2) <= TlowerforColdCC && A2(k,3) >= TupperforColdCC
                sumforcold=sumforcold+A2(k,1)*(TupperforColdCC-TlowerforColdCC)/A2(k,4);
            end
        end
        % calculating area using the formula A = Q/U * 1/LMT
        A = A + 1/LMT * (sumforhot + sumforcold);
        
        [~,TlowerforHotCC] = polyxpoly(xhot,yhot,[temperatureinterval(i) temperatureinterval(i)],[0 1000]);
        [~,TlowerforColdCC] = polyxpoly(xcold,ycold,[temperatureinterval(i) temperatureinterval(i)],[0 1000]);
        [~,TupperforHotCC] = polyxpoly(xhot,yhot,[xhot(1) xhot(1)],[0 1000]);
        [~,TupperforColdCC] = polyxpoly(xcold,ycold,[(xhot(1)-0.00001) (xhot(1)-0.00001)],[0 1000]);
       
        TlowerforHotCC = max(TlowerforHotCC);
        TupperforHotCC = max(TupperforHotCC);
        TlowerforColdCC = max(TlowerforColdCC);
        TupperforColdCC = max(TupperforColdCC);
        
        LMT = ((TupperforHotCC - TupperforColdCC)-(TlowerforHotCC - TlowerforColdCC))...
            /log((TupperforHotCC - TupperforColdCC)/(TlowerforHotCC - TlowerforColdCC));
        
        sumforhot = 0;
        sumforcold = 0;
        
        for j = 1:numel(A1(:,1))
            if A1(j,2) >= TupperforHotCC && A1(j,3) <= TlowerforHotCC
                sumforhot=sumforhot+A1(j,1)*(TupperforHotCC-TlowerforHotCC)/A1(j,4);
            end
        end
        
        for k = 1:numel(A2(:,1))
            if A2(k,2) <= TlowerforColdCC && A2(k,3) >= TupperforColdCC
                sumforcold=sumforcold+A2(k,1)*(TupperforColdCC-TlowerforColdCC)/A2(k,4);
            end
        end
        
        A = A + 1/LMT * (sumforhot + sumforcold);
        
    else
        LMT = ((TupperforHotCC - TupperforColdCC)-(TlowerforHotCC - TlowerforColdCC))...
            /log((TupperforHotCC - TupperforColdCC)/(TlowerforHotCC - TlowerforColdCC));
        
        sumforhot = 0;
        sumforcold = 0;
        for j = 1:numel(A1(:,1))
            if A1(j,2) >= TupperforHotCC && A1(j,3) <= TlowerforHotCC
                sumforhot=sumforhot+A1(j,1)*(TupperforHotCC-TlowerforHotCC)/A1(j,4);
            end
        end
        
        for k = 1:numel(A2(:,1))
            if A2(k,2) <= TlowerforColdCC && A2(k,3) >= TupperforColdCC
                sumforcold=sumforcold+A2(k,1)*(TupperforColdCC-TlowerforColdCC)/A2(k,4);
            end
        end
        
        A = A + 1/LMT * (sumforhot + sumforcold);      
        
    end
end
    number_seg = Wheretosplit(A1,A2);
    
    hold on
    title(['Approximate Heat Exchanger Area Required = ',num2str(round(A,2)), ' m^{2}'...
        newline 'Number of Segments Involved = ', num2str(number_seg)],'FontSize', 23)
end









function [output, output1, output2] = HotCC(A)

% HotCC only works for BCC, it plots the hot composite curve
% input A is a matrix:
%   [CPstream1 Tsupply Ttarget; CPstream2 Ts Tt; CPstream3 Ts Tt; ...]


% sort unique temperatures in descending order
temperatureinterval = sort(unique(A(:,2:3)),'descend');

% get enthalpy change for each temperature interval
EnthalpyChange=[];
for i = 2:numel(temperatureinterval)
    T = temperatureinterval(i);
    CP=[];
    for j = 1:numel(A(:,1))
        if A(j,2) > T && A(j,3) <= T % check if a particular stream is within the T interval
            CP = [CP A(j,1)];
        end 
    end
    EnthalpyChange = [EnthalpyChange sum(CP)*(temperatureinterval(i-1) - T)];
end

enthalpysum = sum(EnthalpyChange);

% get data for plotting
output = [enthalpysum - EnthalpyChange(1)];
for i = 2:numel(EnthalpyChange)-1
    enthalpy = output(i-1)-EnthalpyChange(i);
    output = [output enthalpy];
end

grid on
hold on
hotCC = plot([enthalpysum output 0],temperatureinterval,'Marker','.',...
    'MarkerSize',25,'LineWidth', 2.6, 'DisplayName', 'Hot CC', ...
    'Color',[0.929 0.490 0.192]);

output1 = [enthalpysum output 0];
output2 = temperatureinterval;

legend(hotCC, 'Hot CC')

end




function [output, output1, output2] = ColdCC(A)

% ColdCC only works for BCC, it plots the cold composite curve
% input A is a matrix:
%   [CPstream1 Tsupply Ttarget; CPstream2 Ts Tt; CPstream3 Ts Tt; ...]

% sort unique temperatures in ascending order
temperatureinterval = sort(unique(A(:,2:3)),'ascend');

% get enthalpy change for each temperature interval
EnthalpyChange=[];
for i = 2:numel(temperatureinterval)
    T = temperatureinterval(i);
    CP=[];
    for j = 1:numel(A(:,1))
        if A(j,2) < T && A(j,3) >= T % check if a particular stream is within the T interval
            CP = [CP A(j,1)];
        end  
    end
    EnthalpyChange = [EnthalpyChange sum(CP)*(temperatureinterval(i-1)-T)];
end

% get data for plotting
output = [0 - EnthalpyChange(1)];
for i = 2:numel(EnthalpyChange)-1
    enthalpy = output(i-1)-EnthalpyChange(i);
    output = [output enthalpy];
end

grid on
hold on
xlabel('Heat Flow (kW)', 'FontSize', 20)
ylabel('Temperature (^{o}C)', 'FontSize', 20)

plot([0 output (-1*sum(EnthalpyChange))],temperatureinterval,'Marker',...
    '.','MarkerSize',25,'LineWidth', 2.6, 'DisplayName', 'Cold CC',...
    'Color',[0.267 0.449 0.769]);

output1 = [0 output (-1*sum(EnthalpyChange))];
output2 = temperatureinterval;

set(gca,'FontSize',18)
end



function output = forArea(A1,A2)

% forArea only works for BCC, input A1 for hot streams in the format
%   [CPstream1 Ts Tt; CPstream2 Ts Tt;...]

%input A2 for cold streams in the format
%   [CPstream1 Ts Tt; CPstream2 Ts Tt;...]

output = unique([HotCC(A1),ColdCC(A2)]);

end




function output = Wheretosplit(A1,A2)

% Wheretosplit only works for BCC, it indicates the segments involved in area targeting

% input A1 for hot streams in the format
%   [CPstream1 Ts Tt;CPstream2 Ts Tt;...]

% input A2 for cold streams in the format
%   [CPstream1 Ts Tt;CPstream2 Ts Tt;...]

output = unique([HotCC(A1),ColdCC(A2)]);
Hot = HotCC(A1);
Cold = ColdCC(A2);

A = [Hot Cold];

number_seg = numel(unique(A))+1;

for i=1:numel(output)
    hold on
    segments = plot([output(i) output(i)],[0 230],'LineStyle','--','LineWidth',1,'Color',[0.8 0.1 0.1]);
    set(get(get(segments,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

output = number_seg;

end













