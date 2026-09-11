# AMOC weakening → European crop yields and household energy demand

Tesi, Marco Bova — Università Bocconi.

Stima una *response function* su clima **osservato** (E-OBS, giornaliero), poi la applica
a due futuri indicizzati per indebolimento dell'AMOC (ΔSv):

| ramo | forcing serra | hosing | ΔSv raggiunto |
|---|---|---|---|
| **NAHosMIP** `u03-hos` | nessuno (fondo piControl) | sì, 0.3 Sv | fino a −14.5 Sv |
| **ISIMIP3b** ssp126 / ssp370 | sì | **no** | −3.5 / −5.5 Sv |

ISIMIP **non** è "hosing + CO₂": è uno scenario di emissioni in cui l'AMOC si indebolisce
come risposta normale. Confrontando i due rami **allo stesso ΔSv** si vede quanto della
risposta futura assomiglia a quella dell'AMOC puro.

Outcome: **resa** t/ha (CropStatHarm, NUTS3) e **energia domestica pro capite** (Eurostat,
paese, gas + elettricità).

## Metodo delta

I modelli climatici hanno bias climatologico, quindi il loro output non entra mai
direttamente in una response function empirica. Entra solo la **differenza**:

```
1. modello:  Δ = perturbato − controllo, per cella nativa e per MESE SOLARE
             (additivo in K per la temperatura, rapporto per la precipitazione)
2. regrid:   Δ interpolato bilinearmente sui centri delle celle E-OBS
3. applica:  TG = tg_osservato[cella, giorno] + Δ[cella, mese(giorno)]
4. SOLO ORA: GDD / HDD / CDD / Heat / Frost sul campo GIORNALIERO perturbato
5. impatto:  Δln(outcome) = Σ β · (X_scenario − X_storico)
```

Il passo 3→4 è il punto: gli indicatori sono soglie non lineari (GDD 5/28 °C, dead-band
HDD 15 / CDD 24). Spostare un totale di degree-days già sommato assumerebbe che gli stessi
giorni attraversino le stesse soglie — falso proprio negli scenari che interessano.

**Limite**, vero su entrambi i rami: un delta mensile sposta la media di ogni mese ma non
la **forma** della distribuzione giornaliera (verificato: sd e skewness dentro ogni coppia
cella-mese restano identiche a 1e−15). Più ondate di calore o varianza estiva più larga non
sono rappresentate. ISIMIP è giornaliero e potrebbe evitarlo: la scelta di aggregarlo
comunque a mensile è deliberata, per non far differire i due rami **anche per metodo**.

## Fasi

| fase | cosa | entry point |
|---|---|---|
| 0 | assemblaggio dati grezzi, rese CropStatHarm su geometria NUTS 2016 | `agriculturegeografical.R` |
| 1 | meteo storico da E-OBS + pannelli di stima | `cropweather_eobs_nuts3.R`, `energyweather_eobs_country.R`, `gdd_window_daily.R` |
| 2 | stima delle response function | `regression.R` |
| 3 | campi delta di scenario | `anomaly/`, `amoc_bin_fields.R`, `isimip_bin_fields.R` |
| 4 | scenario replay (delta applicato a E-OBS) | `scenario_engine.R` + 5 runner |
| 5 | impatto Σ β·ΔX | `amoc_impact.R` + 7 script |

`weather_indicators.R` definisce GDD/Heat/Frost/HDD/CDD **una volta sola** ed è usato sia in
stima sia in replay: se il regressore non è costruito identicamente, β non è trasportabile.
`regression.R` è la fonte della specifica — gli script di Fase 5 la leggono con `grab()`,
nessuna formula è mai riscritta.

## ⚠️ Cosa gira da un clone pulito

Il repo contiene gli **script** e gli **output tabellari**, non i dati grezzi (~58 GB) né gli
intermedi rigenerabili.

| | stato |
|---|---|
| **Fase 2** (`regression.R`) | ✅ gira — i pannelli `9.` e `16.` sono in repo |
| **Plot** (`plot_impact_curve_*.R`) | ✅ girano — leggono solo i `.csv` di impatto, tutti in repo |
| **Fase 5** (`amoc_impact.R` e discendenti) | ❌ **non gira** — servono i `scenario_*_bins_*.csv.gz` |
| **Fase 4** (`scenario_replay_*.R`) | ❌ non gira — servono i `.nc` dei campi delta |
| **Fasi 0/1/3** | ❌ non girano — servono i dati grezzi |

`amoc_impact.R` è la base che ogni altro script di Fase 5 fa `source()`, quindi la mancanza
dei `.csv.gz` blocca l'intera fase. Sono esclusi di proposito: la tabella crop è ~4.2 M righe,
oltre il limite per file di GitHub.

Per rigenerare, nell'ordine, avendo i dati grezzi:

```bash
Rscript Amoc/Code/anomaly/nahos_anomaly_cdo\ \(3\).sh   # + i due script EC-Earth3 / IPSL
Rscript Amoc/Code/amoc_bin_fields.R                     # -> amoc_{bin,year}_fields_*.nc
Rscript Amoc/Code/isimip_bin_fields.R                   # -> isimip_{bin,year}_fields_*.nc
Rscript Amoc/Code/scenario_replay_bins.R                # -> scenario_bins_*.csv.gz     (ore)
Rscript Amoc/Code/scenario_replay_isimip_bins.R         # -> scenario_isimip_bins_*.csv.gz
Rscript Amoc/Code/amoc_impact.R                         # -> amoc_impact_*.csv
```

`ISIMIP_SCEN=ssp370` su qualunque script del ramo ISIMIP produce le varianti `_ssp370`;
senza variabile d'ambiente il default è `ssp126` e i nomi restano nudi.

Ogni runner di Fase 4 apre con un **zero-delta self-check**: con Δ=0 e rapporto=1 il motore
deve riprodurre *esattamente* i file storici `8.` / `11.` / `13.`, altrimenti `stop()`. Non è
un warning: ogni numero a valle è una differenza di questi indicatori, quindi un difetto del
motore li distorcerebbe tutti in silenzio.

## ⚠️ Il risultato crop non è un risultato

Il segno dell'effetto AMOC sulle rese dipende interamente da `β_gdd < 0`, che non è
interpretabile. Tre specificazioni, tre fallimenti diversi:

| spec | finestra | varianza identificante | esito |
|---|---|---|---|
| **A** | Mar–Lug fissa | sì | β preciso ma confuso col disallineamento fenologico |
| **C** | tempo termico | **0.0 %** | non identificato per costruzione |
| **D** | C + `n_day` | 31–37 % | identificato, non significativo in 3 colture su 4 |

`crop_dropgdd_bound.R` quantifica l'esposizione: togliendo `gdd_mj` il segno si ribalta in
**37 bin su 40**. Da riportare come *scelta di specificazione*, non come finding. Le
componenti restano informative; il netto no. Il ramo **energy** non ha questo problema —
le sue finestre sono costruzioni di calendario che corrispondono a cicli di consumo reali.

Dettagli in `Amoc/Methodology Explanation/` (Appendici A–J).
