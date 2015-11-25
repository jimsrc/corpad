clear
%
DIRECTORIO_DATA		= '../'
ARCHIVO_INPUT_PLAS 	= sprintf('%sINPUT_PLA.inp', DIRECTORIO_DATA);
ARCHIVO_SALIDA 		= 'sqr_displacements.dat';
RIGIDITY		= 1.0
%
data_inpu_plas 	= load(ARCHIVO_INPUT_PLAS);
N_GYROPERIODS 	= data_inpu_plas(3);
N 		= data_inpu_plas(4);
T_MAX_PARTIAL	= N_GYROPERIODS / N;
%
for i = 1:N
	fprintf('archivo #%03d \n', i-1)
	FILENAME = sprintf('%splas_%4.4fGV_tmax_%03d.out', DIRECTORIO_DATA, RIGIDITY, i-1);
	data = load(FILENAME);

	dx2 = data(:,1).*data(:,1);
	dy2 = data(:,2).*data(:,2);
	dz2 = data(:,3).*data(:,3);
	aux_perp2 = dx2 + dy2;

	dr_perp2(i) = mean(aux_perp2);
	dr_para2(i) = mean(dz2);
	time_max(i) = i * T_MAX_PARTIAL;
end
%---------------encabezado
dlmwrite(ARCHIVO_SALIDA, '# $1: time_max_partial 	[1/omega]', 'delimiter', '')
dlmwrite(ARCHIVO_SALIDA, '# $2: <dr_perp^2> 	[km^2]', 'delimiter', '', '-append')
dlmwrite(ARCHIVO_SALIDA, '# $3: <dr_parall^2> 	[km^2]', 'delimiter', '', '-append')
%---------------guarda data
DISPLACEMENTS = cat(2, time_max', dr_perp2', dr_para2');
dlmwrite(ARCHIVO_SALIDA, DISPLACEMENTS, 'delimiter', ' ', 'precision', 12, '-append')
%
%%
