clc;clear;
id=fopen('PeriodicTable.txt'); fgetl(id);data=textscan(id,'%f %s %s %s %s %f %f %s'); fclose(id);
AtomSymbol=data{1,2}; RadioType=data{1,5};EA=data{1,6};IP=data{1,7};
Electronegativity=(EA+IP)/2;Isuse=data{1,8};
ind= ismember(Isuse,'1') & ismember(RadioType,'Non-Radioactive');
AtomSymbol=AtomSymbol(ind);Electronegativity=Electronegativity(ind);

load mp_index.mat;
load DirectGap.mat;% direct energy gap : 271*1 cell
id=fopen('mp_filter_results.txt','w');aa=[];
for ii=1:length(mp_index)
    fprintf('ii=%d of %d\n',ii,length(mp_index));
    material_id=mp_index(ii).material_id;   icsd_id=mp_index(ii).icsd_ids; band_gap=mp_index(ii).band_gap;
    e_above_hull=mp_index(ii).e_above_hull; formula=mp_index(ii).formula;  atomType=formula(1:2:end);
    spacegroup=mp_index(ii).spacegroup;
    [isMem,atomNum]=ismember(atomType,AtomSymbol);
    
    A = sum(isMem)==length(atomType);% A: formula in allowed elements
    B = ~isempty(icsd_id);% B: have ICSD
    C = e_above_hull<=0.02;% C: Ehull<=0.02
    if A
        K=Electronegativity(atomNum); N=cell2mat(formula(2:2:end))';
        temp=prod(K.^N);Kcomp=temp^(1/sum(N));
        D = Kcomp>=4 && Kcomp<=6;% D: 4<=Kcomp<=6
    else
        D = false;
    end
    E = band_gap>0 && band_gap<=3;% E: 0<Egap<=3 eV
    F = ismember(material_id,direct);% F: have direct band gap
    
    if A && B && C && D && E && F
        fa=mp_index(ii).formula;
        spacegroup=mp_index(ii).spacegroup{4};gap=mp_index(ii).band_gap;
        eah=mp_index(ii).e_above_hull;vol=mp_index(ii).volume;
        fprintf(id,'%s ',material_id);
        f1=formula(1:2:end);f2=formula(2:2:end);
        for jj=1:length(f1)
            fprintf(id,'%s',f1{jj});
            if f2{jj}>1
                fprintf(id,'%d',f2{jj});
            end
        end
        fprintf(id,' %s ',spacegroup);
        fprintf(id,'%.4f ',eah);
        fprintf(id,'%.4f ',gap);
        fprintf(id,'%.4f ',Kcomp);
        ids=int2str(mp_index(ii).icsd_ids{1});icsd='icsd_000000';icsd(end-length(ids)+1:end)=ids;
        file=['icsd_poscar/',icsd,'.vasp'];
        ge=import_poscar(file);Nsites=sum(ge.atomcount);
        fprintf(id,'%.0f ',Nsites);
        fprintf(id,'%.4f ',vol);
        if rem(Nsites/sum(cell2mat(f2)),1)~=0 % error('Wrong Nsites !!!');
            fprintf(id,'%s','Nsites/(Nformula) ~= integer, please be aware !!!!!!!!!!!!!!!!!');
        end
        fprintf(id,'\n');
    end
  
end
fclose all;