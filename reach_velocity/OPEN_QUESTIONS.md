# Open questions for the field team

Things the analysis needs from the people who ran the surveys and the instruments. Each one says what was found in the data, why it matters, and what would settle it. None of them changes the Middle-reach results; the first affects one figure panel, the rest affect the other two reaches.

## 1. The ADCP: what does 'Vb' measure?

The export (`data_raw/25_ADCPData.csv`) has one field, `sontek_Vb`, in 0.01 m/s steps. The nearest surveyed profile (1.7 m from the instrument) has its bed 4 cm below the stage sensor's datum, so the instrument sits under about 0.8 m of water at 0.75 m of stage. At that stage it reads 0.075 m/s, while the profile beside it gives a column mean of 0.24 m/s; 0.075 m/s is what that profile's log law gives 3 cm above the bed. Its stage response is on the ridge the transects give (exponent 1.1 to 1.3 on stage − 0.45 m), so the shape check stands either way. The level check (Figure 3b of the paper) assumes Vb is one slanted beam's along-beam velocity, which for a 25° beam is 0.42 of the along-stream speed; that assumption explains the level to within a tenth (ratio 0.91, an implied beam angle of 22.6°).

Needed: the instrument's model (IQ, Argonaut-SW, other), what `Vb` is (a single beam, a bottom cell, the profile mean), the cell size and blanking distance, and whether it logged its own vertical-beam depth. If Vb is not a single beam, the ADCP goes back to a shape check with the level fitted, and the paper still stands.

## 2. The velocity meter

The make and model of the hand-held meter, and its stated accuracy near 0.02 m/s, for the methods section.

## 3. The 0.2 m/s convention

The source of the 0.2 m/s threshold used in the original notebooks. It now enters only as a comparison (the nest-learned rule replaces it), but it should be cited.

## 4. The 'BT' transects (other reaches)

`0709-BT01`, `0730-BT01`, `0813-BT01` (Lower) and `0708-BT01`, `0729-BT01` (Middle), July and August 2026, have one reading per point at 0.4 of the depth (like the 2025 'Other' rows) and RTK elevations about 7 m above the reach's water surface (595 to 600 ft against 571 to 577 ft for the rest of the Lower). They are set aside. Was this a different survey (bank transect, boat), and was the GNSS set up differently that day? The two Middle ones never entered the analysis (one has no coordinates, one no reading heights).

## 5. RTK sessions without a fix (other reaches)

Across a transect the water surface is flat, so RTK bed plus rod depth should agree between points to a few centimetres. It does in June (1 to 9 cm). On 10, 12 and 19 August 2026 (Upper) the points disagree by 0.4 to 1.2 m, and on 6 August (Lower) by 0.15 to 0.3 m: a float or autonomous solution. The bed shapes are still usable (bed = the transect's water surface minus the rod depth), but those days cannot be tied to the other days' water levels. Was the rover in float, or the base set up differently?

## 6. The Lower reach's water level and the dam

The Lower sensor (ARB112) sits 258 m above the first section, below Flat Rock dam. The 15 June 2026 survey was made during a release pulse (the sensor rose from 0.41 to 1.11 m that morning), which is why that day's depths are 60 cm above the neighbouring days at the same sections. A release schedule for June 2026, if one exists, would say whether the pulse was the dam or a storm. The sensor's record also carries a daytime artefact of 2 to 9 cm (four times the Middle's) and spikes; it needs the Middle's cleaning before use (logbook entry 24).

## 7. The Upper sensor's site

The Upper sensor (ARB017) sits on a different branch of the river network from the surveyed centreline, 463 m from it and 1.4 km from the nearest section, on a damped, lake-like site. On the tied survey days the sections' water level moves about 3.5 m per metre of sensor stage. Is this the intended site for the Upper reach, or is there a sensor closer to the sections?
