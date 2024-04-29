# Notes

## Research

## Queries

#### Gaia sources, SDSS cross-match and photometry
```
SELECT gs.source_id, sdss.clean_sdssdr13_oid, sdss.original_ext_source_id, sdss.angular_distance,
       gs.ra AS gaia_ra, gs.dec AS gaia_dec, gs.phot_g_mean_mag, gs.bp_g, gs.g_rp
FROM gaiadr3.xp_summary xp
JOIN gaiadr3.sdssdr13_best_neighbour sdss ON sdss.source_id = xp.source_id
JOIN gaiadr3.gaia_source gs ON gs.source_id = xp.source_id
WHERE xp.bp_n_contaminated_transits = 0 AND xp.rp_n_contaminated_transits = 0
AND gs.phot_g_mean_mag < 20.5
AND gs.ra > 60 AND gs.ra < 300 AND dec > -1.67 AND dec < 1.67
```
- the last line is a ra/dec filter for Stripe 82

#### SDSS sources - redshift, photometry, specObjId
```
SELECT objId, specObjId, z, zErr, ra AS specRa, dec AS specRa, photoRa, photoDec,
       psfMag_u, psfMag_g, psfMag_r, psfMag_i, psfMag_z,
       psfMagErr_u, psfMagErr_g, psfMagErr_r, psfMagErr_i, psfMagErr_z
FROM SpecPhoto
WHERE psfMag_r < 20
AND ra > 60 AND ra < 300 AND dec > -1.67 AND dec < 1.67
```
- the last line is a ra/dec filter for Stripe 82


## SDSS spectra (full sky) - G < 17
- 146745 total (speclist)
- legacy: 100380
- stellar cluster plates: 3369
- SEGUE-2: 0
- eBOSS: 15703
