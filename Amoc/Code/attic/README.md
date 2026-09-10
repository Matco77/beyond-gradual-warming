# attic

Script superati, tenuti per traccia. **Non fanno parte della pipeline**: nessuno
li sorgenta, nessuno legge i loro output.

| file | perché è qui |
|---|---|
| `Dataset_Preparation.R` | notebook esplorativo iniziale. Legge un CSV FAOSTAT che non elabora mai, definisce `inspect_nc_table()` per ispezionare i NetCDF NAHosMIP, e apre alcuni raster HadGEM3-GC31-LL. Non gira comunque: chiama `inspect_nc_table(g01_tas_1)` alla riga 55, ma `g01_tas_1` è definito alla 74. I suoi due output (`g01_tas_dimensions_table.csv`, `g01_tas_variables_table.csv`) non sono letti da nessuno script. Il file rimanda esso stesso a `agriculturegeografical.R`, che è la versione buona di quel passaggio. |
