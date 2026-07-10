from pathlib import Path

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import ConnectionPatch

PROJECT_ROOT = Path(__file__).resolve().parent
OUT_DIR = PROJECT_ROOT / "figures" / "site_description"
OUT_DIR.mkdir(parents=True, exist_ok=True)

lat_point = 49.712
lon_point = -2.010

fig = plt.figure(figsize=(10, 10)) 
main_ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

main_ax.set_extent([-5, 10, 41, 52], crs=ccrs.PlateCarree())
main_ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#f2f2f2')
main_ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#b3d9ff')
main_ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
main_ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', alpha=0.5)

main_ax.plot(lon_point, lat_point, 'ro', markersize=8, markeredgecolor='white', markeredgewidth=2, transform=ccrs.PlateCarree())

zoom_radius = 0.4 
zoom_extent = [lon_point - zoom_radius, lon_point + zoom_radius, 
               lat_point - zoom_radius, lat_point + zoom_radius]

zoom_ax = fig.add_axes([0.58, 0.58, 0.35, 0.35], projection=ccrs.PlateCarree())
zoom_ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())

zoom_ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#f2f2f2')
zoom_ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#b3d9ff')
zoom_ax.add_feature(cfeature.COASTLINE.with_scale('10m'))
zoom_ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', alpha=0.5)

zoom_ax.plot(lon_point, lat_point, 'ro', markersize=8, markeredgecolor='white', markeredgewidth=2, transform=ccrs.PlateCarree())

zoom_ax.spines['geo'].set_edgecolor('red')
zoom_ax.spines['geo'].set_linewidth(2)

con = ConnectionPatch(xyA=(lon_point, lat_point), coordsA=main_ax.transData,
                      xyB=(lon_point - zoom_radius, lat_point), coordsB=zoom_ax.transData,
                      arrowstyle="-", 
                      linestyle="--", 
                      color="#333333",  
                      linewidth=2.5,   
                      alpha=1.0)
fig.add_artist(con)

plt.title("Coordinates :    ", fontsize=22.5)
plt.savefig(OUT_DIR / "alderney_race_map.png", dpi=300, bbox_inches="tight")
plt.show()
