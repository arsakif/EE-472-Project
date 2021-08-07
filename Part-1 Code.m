function [Y_bus, t_solution, number_of_non_zero_elements, V_bus] = e194908_arslan_PF(argin1)
% EE 472 Power System Analysis II Fall 2019
% Term Project Part I: Submitting Bus Admittance Matrix.
% Instructor         : Assist. Prof. Murat GOL
% Course Assistant.  : Mustafa Erdem Sezgin
%===========================================================================
% Akif ARSLAN 1949080
%===========================================================================
% FORM Y_BUS
tic
ieee300 = fopen(argin1,'r');
% Find where bus data starts.
while true
    line_data = fgetl(ieee300);
    if length(line_data) >= 3
        if strcmpi(line_data(1:3),'BUS')
            break
        end
    end
end
% allocate indexing matrix for faster iteration
indexing_matrix = uint16(zeros(1000,1));
shunt_G_and_B = zeros(1000,1);
i = uint16(1);
while true
    line_data = fgetl(ieee300);
    % break the loop when BUS DATA finishes
    if length(line_data) < 50
        % cut unused parts
        indexing_matrix = indexing_matrix(1:i-1);
        shunt_G_and_B = shunt_G_and_B(1:i-1);
        number_of_buses = i-1;
        break
    end
    indexing_matrix(i) = str2double(line_data(1:4));
    shunt_G_and_B(i) = str2double(line_data(107:114)) + str2double(line_data(115:122))*1i;
    i = i + 1;
end
% Go where BRANCH DATA STARTS
while true
    line_data = fgetl(ieee300);
    if length(line_data) >= 6
        if strcmpi(line_data(1:6),'BRANCH')
            break
        end
    end
end
% allocate memory for Y_BUS matrix for faster iteration
Y_BUS = zeros(number_of_buses*4,4);
%Create a 2x2 temporary pair admittance matrix
Y_BUS_temp = zeros(2,2);
k=uint16(1);
while true
    line_data = fgetl(ieee300);
    % break the loop when BRANCH DATA finishes
    if length(line_data) < 50
        if sum(Y_BUS(:,1) == 0,1)
            % Cut unused parts of the Y_BUS
            zero_cut = find(Y_BUS(:,1)==0, 1, 'first');
            Y_BUS = Y_BUS(1:zero_cut-1,:);
        end
        break
    end
    %Use indexing as given in BUS data order.
    Yi = find(indexing_matrix == str2double(line_data(1:4))); % "from" bus
    Yj = find(indexing_matrix == str2double(line_data(6:9))); % "to" bus
    %Get Resistance and Reactance data and turn into line admittance.
    %R = str2double(line_data(20:29));
    %X = str2double(line_data(30:40));
    line_admittance = 1/(str2double(line_data(20:29)) + str2double(line_data(30:40))*1i);
    %Get line charing B data and divide by 2
    line_charging = (str2double(line_data(41:50))/2)*1i;
    % If there is any tap or phase shifter include their effects.
    switch str2double(line_data(19))
        case 0 % 0 ==> A line.
            Y_BUS_temp(1,1) = line_admittance + line_charging + shunt_G_and_B(Yi);% Yii
            Y_BUS_temp(1,2) = -line_admittance; % Yij
            Y_BUS_temp(2,1) = -line_admittance; %Yji
            Y_BUS_temp(2,2) = line_admittance + line_charging + shunt_G_and_B(Yj); % Yjj
        case {1,2,3} % 1,2,3 ==> There is a tap changer.
            tap = str2double(line_data(77:82));
            Y_BUS_temp(1,1) = (line_admittance/(tap^2)) + shunt_G_and_B(Yi);% Yii
            Y_BUS_temp(1,2) = -line_admittance/tap; % Yij
            Y_BUS_temp(2,1) = -line_admittance/tap; % Yji
            Y_BUS_temp(2,2) = line_admittance + shunt_G_and_B(Yj); % Yjj
        case 4 % 4 ==> There is a phase shifter.
            tap = str2double(line_data(77:82));
            p_shift = str2double(line_data(84:90));
            p_shift = p_shift*pi/180;
            Y_BUS_temp(1,1) = line_admittance/tap^2 + shunt_G_and_B(Yi);% Yii
            Y_BUS_temp(1,2) = -line_admittance/(cos(p_shift) - sin(p_shift)*1i); % Conjugate
            Y_BUS_temp(2,1) = -line_admittance/(cos(p_shift) + sin(p_shift)*1i);
            Y_BUS_temp(2,2) = line_admittance + shunt_G_and_B(Yj);
        otherwise
            disp('line information has not found')
    end
    % Once shunt values are used remove them to avoid adding again at next
    % iterations.
    shunt_G_and_B(Yi) = 0;
    shunt_G_and_B(Yj) = 0;
    
    is_Yi_used = logical(sum(any(Y_BUS(:,1:2) == Yi)));
    is_Yj_used = logical(sum(any(Y_BUS(:,1:2) == Yj)));
    % If both Yi and Yj busses are used at previous iteration, don't create
    % new Yii and Yjj. Only create Yij and Yji.
    if is_Yi_used && is_Yj_used
        index_temp = find(sum(Y_BUS(:,[1,2])==Yi,2)==2);
        Y_BUS(index_temp,3) = Y_BUS(index_temp,3) + real(Y_BUS_temp(1,1));
        Y_BUS(index_temp,4) = Y_BUS(index_temp,4) + imag(Y_BUS_temp(1,1));
        
        index_temp = find(sum(Y_BUS(:,[1,2])==Yj,2)==2);
        Y_BUS(index_temp,3) = Y_BUS(index_temp,3) + real(Y_BUS_temp(2,2));
        Y_BUS(index_temp,4) = Y_BUS(index_temp,4) + imag(Y_BUS_temp(2,2));
    else
        % If one of the buses are used to calculate Y_BUS in previous
        % steps, don't creat new entry for it. Add its new admittance value
        % to previous one.
        if is_Yi_used || is_Yj_used
            % Check, which one of Yi and Yj are used before
            if is_Yi_used
                index_temp = find(sum(Y_BUS(:,[1,2])==Yi,2)==2);
                Y_BUS(index_temp,3) = Y_BUS(index_temp,3) + real(Y_BUS_temp(1,1));
                Y_BUS(index_temp,4) = Y_BUS(index_temp,4) + imag(Y_BUS_temp(1,1));
                
                Y_BUS(k,1) = Yj;
                Y_BUS(k,2) = Yj;
                Y_BUS(k,3) = real(Y_BUS_temp(2,2));
                Y_BUS(k,4) = imag(Y_BUS_temp(2,2));
                k = k + 1;
            else % if is_Yj_used
                Y_BUS(k,1) = Yi;
                Y_BUS(k,2) = Yi;
                Y_BUS(k,3) = real(Y_BUS_temp(1,1));
                Y_BUS(k,4) = imag(Y_BUS_temp(1,1));
                k = k + 1;
                
                index_temp = find(sum(Y_BUS(:,[1,2])==Yj,2)==2);
                Y_BUS(index_temp,3) = Y_BUS(index_temp,3) + real(Y_BUS_temp(2,2));
                Y_BUS(index_temp,4) = Y_BUS(index_temp,4) + imag(Y_BUS_temp(2,2));
            end
            % If neighter Yi or Yj bus number is used previously,
            %construct new Yii and Yjj.
        else
            Y_BUS(k,1) = Yi;
            Y_BUS(k,2) = Yi;
            Y_BUS(k,3) = real(Y_BUS_temp(1,1));
            Y_BUS(k,4) = imag(Y_BUS_temp(1,1));
            k = k + 1;
            
            Y_BUS(k,1) = Yj;
            Y_BUS(k,2) = Yj;
            Y_BUS(k,3) = real(Y_BUS_temp(2,2));
            Y_BUS(k,4) = imag(Y_BUS_temp(2,2));
            k = k + 1;
        end
    end
    % Cunstruc Yij = Yji element for any of the cases.
    Y_BUS(k,1) = Yi;
    Y_BUS(k,2) = Yj;
    Y_BUS(k,3) = real(Y_BUS_temp(1,2));
    Y_BUS(k,4) = imag(Y_BUS_temp(1,2));
    k = k + 1;
    
    Y_BUS(k,1) = Yj;
    Y_BUS(k,2) = Yi;
    Y_BUS(k,3) = real(Y_BUS_temp(2,1));
    Y_BUS(k,4) = imag(Y_BUS_temp(2,1));
    k = k + 1;
end
% Close the text file after iteration.
fclose(ieee300);
% Use "sparse" fuction to obtain Y_BUS as a sparse matrix. sparse(i,j,v)
Y_bus = sparse(Y_BUS(:,1),Y_BUS(:,2), Y_BUS(:,3) + Y_BUS(:,4)*1i);
number_of_non_zero_elements = nnz(Y_bus);
% Create a random I vector.
I = rand(number_of_buses,1);
% Calculate V_bus
V_bus = Y_bus\I;
t_solution = toc;

