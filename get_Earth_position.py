from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth

# Earth position in Skycoord
earth_position = get_earth('2021-01-17 12:00:00')

# transform to HCI coordinate
earth_HCI_coord = earth_position.transform_to(frames.HeliocentricInertial)
print('Earth in HCI coordinate: ', earth_HCI_coord)
