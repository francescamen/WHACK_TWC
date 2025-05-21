close all
clc
clear

%% configuration
BW = 80;  % bandwidth

routers_csi = string([3]);
routers_csi_num = [3];

routers = (1:length(routers_csi));

folder = "";
VSA='0003';
VSB='000';

subfolder_name = "020822";

d = dir(strcat("../../traces/", subfolder_name));
mkdir(strcat("../../mat_files/", subfolder_name));

d_cell = struct2cell(d);
name_folder = string(d_cell(1,3:end).')

for id_point = 1:length(name_folder)
    mkdir(strcat("../../mat_files/", subfolder_name, "/", name_folder(id_point)))

    for id_router = routers
        
        [id_point, id_router]
        
        VSB(4) = lower(dec2hex(routers_csi_num(id_router)));
        VSB
        
        FILEA = strcat("../../traces/", subfolder_name, "/", ...
            name_folder(id_point), "/trace", routers_csi(id_router), ...
            ".pcap")
        [~, ~, ~, ~, ~, ~, ~, csi_store, toa_packets ,~] = load80MHZ_no_correz(FILEA, BW, VSA);
        packets = length(csi_store);
        K = length(csi_store{1,1}.core{1,1}.nss{1,1}.data);
        N = csi_store{1,1}.mask_to_process;
        
        csi_data = zeros(packets, K, N);
        
        for ii = 1:packets
            for jj = 1:N
                csi_data(ii,:,jj) = csi_store{1,ii}.core{1,jj}.nss{1,1}.data;
            end
        end
        save(strcat("../../mat_files/", subfolder_name, "/", ...
            name_folder(id_point), "/csi_data"), "csi_data", "toa_packets")
    end
end
