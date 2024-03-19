import json
from collections import defaultdict

import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from geopy.geocoders import Nominatim
from matplotlib.cm import Dark2
from matplotlib.cm import Set1
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.ops import nearest_points
from shapely.prepared import prep
from tqdm import tqdm
from h3 import h3
import pyproj
from functools import partial
from shapely.ops import transform
from pyproj import Transformer
import matplotlib.pyplot as plt
import contextily as ctx
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def get_city_boundary(place_name, city_name='Jaipur District', boundary_name='Jaipur Municipal Corporation', plot=True):
    """
    Retrieves the administrative boundary of a specified city and optionally plots it.

    Parameters:
    - place_name (str): The name of the city for which to retrieve the boundary, including country for disambiguation.
    - city_name (str): The specific administrative name of the city district to filter the search. Default is 'Jaipur District'.
    - boundary_name (str): The name of the boundary level to retrieve, typically representing a municipal corporation or a similar administrative level. Default is 'Jaipur Municipal Corporation'.
    - plot (bool): If True, the function plots the boundary. Default is True.

    Returns:
    - city_geometry (shapely.geometry): The unified geometry object representing the city's administrative boundary.
    """
    # Get shapefile of the city based on administrative boundary
    city = ox.geometries_from_address(place_name, tags={'boundary': 'administrative'})
    city = city[city.name == city_name]

    # Get the city boundary at the specified boundary level
    city_boundary = ox.geometries_from_polygon(city.geometry.unary_union, tags={'boundary': 'administrative'})
    city_boundary = city_boundary[city_boundary.name == boundary_name]

    # Plot the city boundary if plot is True
    if plot:
        city_boundary.plot()

    # Get the unified geometry of the city boundary
    city_shape = city_boundary.geometry.unary_union
    city_geometry = gpd.GeoSeries([city_shape], crs="EPSG:4326")

    return city_shape, city_geometry


def reset_crs(gdf, crs='EPSG:4326'):
    """Resets CRS EPSG:4326 for GeoPandas DataFrame"""
    gdf = gdf.to_crs(crs)
    return gdf


def set_local_crs(gdf, crs='EPSG:32643'):
    """Sets local CRS EPSG:32643 for GeoPandas DataFrame"""
    gdf = gdf.to_crs(crs)
    return gdf


def create_grid(geometry, grid_cell_size=1000, plot=False, native_crs='EPSG:32643'):
    """
    Partitions a geometry into a grid and returns a GeoDataFrame.

    Parameters:
    - geometry: GeoSeries or GeoDataFrame containing the geometry.
    - grid_cell_size: int, the length/width of the grid cell in meters.
    - plot_grid: bool, whether to plot the grid or not.
    - native_crs: str, the coordinate reference system for the grid.

    Returns:
    - grid: GeoDataFrame, the geometry partitioned into a grid with crs 'EPSG:4326'.
    """

    def grid_bounds(geom, delta):
        minx, miny, maxx, maxy = geom.bounds
        gx, gy = np.arange(start=minx, stop=maxx, step=delta), np.arange(start=miny, stop=maxy, step=delta)
        grid = []
        for i in range(len(gx) - 1):
            for j in range(len(gy) - 1):
                poly_ij = Polygon([[gx[i], gy[j]], [gx[i], gy[j + 1]], [gx[i + 1], gy[j + 1]], [gx[i + 1], gy[j]]])
                grid.append(poly_ij)
        return grid

    # Method to grid a shape file
    def partition(geom, delta):
        prepared_geom = prep(geom)
        grid = list(filter(prepared_geom.intersects, grid_bounds(geom, delta)))
        return grid

    # Ensure city_geometry is a GeoSeries
    if isinstance(geometry, gpd.GeoDataFrame):
        geom = geometry.geometry.iloc[0]
    else:
        geom = geometry.iloc[0]

    # Partition the geometry into a grid
    grid = partition(geom, grid_cell_size)
    grid = gpd.GeoSeries(grid)

    if plot:
        fig, ax = plt.subplots(figsize=(10, 10))
        grid.boundary.plot(ax=ax)
        gpd.GeoSeries([geom]).boundary.plot(ax=ax, color="red")
        plt.show()

    # Converting geo series grid to a GeoDataFrame with additional columns
    df = pd.DataFrame(grid, columns=['geometry']).reset_index(drop=True)
    df['geometry'] = df['geometry'].astype(str)
    df['geometry'] = df['geometry'].map(lambda x: x.lstrip('POLYGON (').rstrip(')'))
    df['POI'] = gpd.GeoSeries(grid.centroid, crs=native_crs).to_crs('EPSG:4326')
    grid = gpd.GeoDataFrame(df, geometry=polygons_from_geo_series(df['geometry']), crs=native_crs)
    grid = reset_crs(grid)

    return grid

def create_hex_grid(geometry, plot=False, native_crs='EPSG:32643'):

    """
    Partitions a geometry into H3 hexagons, calculates centroids in local CRS,
    and returns a GeoDataFrame with geometries and centroids in 'EPSG:4326'.

    Parameters:
    - geometry: GeoSeries containing the geometry in native_crs.
    - plot: bool, whether to plot the grid or not.
    - native_crs: str, the original coordinate reference system for the geometry.
    - local_crs: str, the local coordinate reference system for calculating centroids.

    Returns:
    - hex_grid: GeoDataFrame, the geometry partitioned into H3 hexagons with centroids, both in 'EPSG:4326'.
    """
    resolution = 8  # Approx. 1 km^2 per hexagon

    # Transform to WGS 84 for H3 processing
    geometry_wgs84 = geometry.to_crs(epsg=4326)

    # Create Transformers for CRS transformations
    transformer_to_local = Transformer.from_crs("EPSG:4326", native_crs, always_xy=True)
    transformer_to_wgs84 = Transformer.from_crs(native_crs, "EPSG:4326", always_xy=True)

    # Function to perform transformation using the transformers
    def transform_geom(geom, transformer):
        return transform(lambda x, y: transformer.transform(x, y), geom)

    # Create an empty list to store hexagon polygons and their centroids
    hex_polygons = []
    centroid_points = []

    for geom in geometry_wgs84:
        if not geom.is_valid:
            continue

        # Get hexagons covering the geometry
        hexes = h3.polyfill(geom.__geo_interface__, resolution, geo_json_conformant=True)
        for hex in hexes:
            # Convert each hex to a Polygon
            polygon = Polygon(h3.h3_to_geo_boundary(hex, geo_json=True))
            hex_polygons.append(polygon)

            # Project hex to local CRS, calculate centroid, then project back to WGS 84
            local_polygon = transform_geom(polygon, transformer_to_local)
            centroid_local = local_polygon.centroid
            centroid_wgs84 = transform_geom(centroid_local, transformer_to_wgs84)

            centroid_points.append(centroid_wgs84)

    # Create a GeoDataFrame from the polygons and centroids
    hex_grid = gpd.GeoDataFrame(geometry=gpd.GeoSeries(hex_polygons, crs='EPSG:4326'))
    hex_grid['POI'] = gpd.GeoSeries(centroid_points, crs='EPSG:4326')

    if plot:
        ax = hex_grid.plot(color='blue', alpha=0.3, edgecolor='k')
        geometry.to_crs(epsg=4326).plot(ax=ax, color='red', alpha=0.5)
        hex_grid['POI'].plot(ax=ax, color='yellow', marker='o', markersize=5)
        plt.axis('off')
        plt.show()

    return hex_grid

def get_address_from_coordinates(latitude, longitude):
    geolocator = Nominatim(user_agent="geoapiExercises")
    try:
        location = geolocator.reverse((latitude, longitude), exactly_one=True)
        return location.address
    except Exception as e:
        print(f"Error occurred: {e}")
        return None


def polygons_from_geo_series(geo_series):
    """"Returns a GeoPandas DataFrame with Shapely Polygons from GeoSeries"""

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def xy_list_from_string(s):
        # 'x y x y ...' -> [[x, y], [x, y], ...]
        return list(chunks([float(i.strip(',')) for i in s.split()], 2))

    def poly(s):
        """ returns shapely polygon from point list"""
        ps = xy_list_from_string(s)
        return Polygon([[p[0], p[1]] for p in ps])

    polygons = [poly(r) for r in geo_series]

    return polygons


def get_nearest_nodes(gdf, graph):
    """Returns the nearest node for each origin and destination"""
    gdf['node'] = ''
    for index, row in tqdm(gdf.iterrows(), total=len(gdf)):
        gdf.loc[index, 'node'] = ox.distance.nearest_nodes(graph, X=row.geometry.centroid.x, Y=row.geometry.centroid.y)
    return gdf


def get_walking_route(graph, trips):
    """"Find the shortest path between source and destination node on street network"""
    trips['walking_distance'] = ''
    for index, poi in tqdm(trips.iterrows(), total=len(trips)):
        # find the shortest path between nodes
        if poi.destination:
            route = ox.shortest_path(graph, poi.node, poi.destination, weight="travel_time")
            # route length in meters
            edge_lengths = ox.utils_graph.get_route_edge_attributes(graph, route, "length")
            trips.loc[index, 'walking_distance'] = round(sum(edge_lengths))
        else:
            pass
    return trips


def nearest_sap(origins, destinations):
    """Returns a dataframe with the nearest SAP"""

    def near(point, pts=None):
        # find the nearest point and return the corresponding Place value
        nearest = destinations.geometry == nearest_points(point, pts)[1]
        return destinations.loc[nearest, ['geometry', 'node', 'sap_name', 'routes']]

    destinations.rename(columns={'name': 'sap_name'}, inplace=True)
    destinations.rename(columns={'name': 'sap_name'}, inplace=True)
    pts = destinations.geometry.unary_union
    origins.reset_index(inplace=True, drop=True)
    # List of new columns to be added
    columns_to_add = ['nearest_sap', 'destination', 'sap_name', 'routes']
    # Adding each column with an empty string as default value
    for column in columns_to_add:
        origins[column] = ''
    # Iterate through origins i.e. POIs (grid centroid)
    for index, row in tqdm(origins.iterrows(), total=len(origins)):
        # Get the nearest point and its linked node
        nearest_sap = near(row.geometry.centroid, pts)
        # select the nearest point if two points are selected, keep the first one (nearest)
        if len(nearest_sap) > 1:
            nearest_sap = nearest_sap.loc[nearest_sap.index[0], :]
        elif len(nearest_sap) == 1:
            origins.iloc[index, 4] = nearest_sap.geometry.item()
            origins.iloc[index, 5] = nearest_sap.node.item()
            origins.iloc[index, 6] = nearest_sap.sap_name.item()
            origins.iloc[index, 7] = nearest_sap.routes.item()

    return origins


def plot_choropleth_sap_poi(gdf, value_column, sap, color_scale="Viridis", title="Choropleth Map"):
    """
    Plots a choropleth map using Plotly from a GeoDataFrame with polygon geometries.

    Parameters:
    - gdf: GeoDataFrame with polygon geometries.
    - value_column: string, the name of the column in gdf to use for coloring the choropleth.
    - color_scale: string, the color scale to use for the choropleth. Default is 'Viridis'.
    - title: string, the title of the map.
    """
    # Ensure the GeoDataFrame is in EPSG:4326 for compatibility with web mapping
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs(epsg=4326)

    # Convert the GeoDataFrame to GeoJSON for Plotly
    geojson = json.loads(gdf.astype(str).to_json())

    # Plotting the choropleth map
    fig = px.choropleth(gdf,
                        geojson=gdf.geometry.__geo_interface__,
                        locations=gdf.index,
                        color=value_column,
                        color_continuous_scale=color_scale,
                        hover_data=['address'],
                        projection="mercator", )

    if sap.crs != 'EPSG:4326':
        sap = sap.to_crs(epsg=4326)
    sap['lon'] = sap.geometry.x
    sap['lat'] = sap.geometry.y
    # Adding scatter plot on the same figure with hover information
    scatter = go.Scattergeo(
        lon=sap['lon'],
        lat=sap['lat'],
        text=sap['sap_name'],  # This line adds hover text from the 'hover_info' column
        mode='markers',
        marker=dict(
            size=8,
            color='red',
            line=dict(width=1, color='rgba(102, 102, 102)'),
        ),
        name='Service Access Points'
    )
    fig.add_trace(scatter)
    fig.update_geos(fitbounds="locations")
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0}, title=title)
    fig.show()


def assign_route_numbers_based_on_connectivity(G, all_stops_gdf):
    """
    Assigns a route number to bus stops based on their presence on the street network
    between any two stops (including terminal stops) on the route.

    Parameters:
    - G: The street network graph obtained via osmnx.
    - all_stops_gdf: GeoDataFrame of all bus stops with 'sap_name', 'node', 'geometry'.
    - route_name: The name (or number) of the route.
    - routes_dict: Dictionary mapping route names to lists of 'sap_name's.

    Returns:
    - Updated routes_dict with 'sap_name's added for stops between important stops.
    """

    def subset_route_df(all_stops_gdf, route_dict, route_name):
        """Returns a dataframe of bus stops for given route name based on sap_name"""
        route_df = all_stops_gdf[all_stops_gdf['sap_name'].isin(route_dict[route_name])].reset_index(drop=True)
        return route_df

    routes_dict = get_routes_dict(bus_stop=all_stops_gdf)

    route_list = list(routes_dict.keys())
    route_list = [ele for ele in route_list if ele != 'no_route']
    for route_name in route_list:

        df_stops = subset_route_df(all_stops_gdf, routes_dict, route_name)

        # Ensure the all_stops_gdf is in the correct CRS
        if all_stops_gdf.crs.to_string() != G.graph['crs']:
            all_stops_gdf = all_stops_gdf.to_crs(G.graph['crs'])

        # Extract node IDs for important stops
        important_stop_ids = df_stops['node'].tolist()

        # Initialize counter for allocated stops
        counter = 0

        # Iterate through pairs of important stops
        for i in range(len(important_stop_ids) - 1):
            start_stop = important_stop_ids[i]
            end_stop = important_stop_ids[i + 1]

            # Compute the shortest path between start_stop and end_stop
            try:
                shortest_path = nx.shortest_path(G, start_stop, end_stop, weight='length')
                # Map node IDs to sap_name and update routes_dict
                for node in shortest_path:
                    if node in all_stops_gdf['node'].values:
                        sap_name = all_stops_gdf.loc[all_stops_gdf['node'] == node, 'sap_name'].values[0]
                        if sap_name not in routes_dict[route_name]:
                            routes_dict[route_name].append(sap_name)
                            if sap_name in routes_dict['no_route']:
                                routes_dict['no_route'].remove(sap_name)
                            counter += 1
            except nx.NetworkXNoPath:
                print(f"No path found between {start_stop} and {end_stop} for route {route_name}")

        print(f'{counter} stops allocated to route {route_name}')
    return routes_dict


def update_destinations_with_routes(destinations, routes_dict):
    """
    Updates the 'routes' column in the destinations DataFrame based on matching 'sap_name' with the routes_dict.

    Parameters:
    - destinations: DataFrame containing columns 'sap_name' and 'routes'.
    - routes_dict: Dictionary with route names as keys and lists of sap_names as values.
    """
    # Iterate through the destinations DataFrame
    for index, row in destinations.iterrows():
        # Initialize an empty list to store routes for the current stop
        routes_for_stop = []

        # Check each route in routes_dict to see if the current stop is included
        for route_name, stops in routes_dict.items():
            if row['sap_name'] in stops:
                # If the stop is part of the route, add the route name to the list
                routes_for_stop.append(route_name)

        # Update the 'routes' column with a comma-separated list of routes for the current stop
        # Only update if there are any routes found; otherwise, leave as None or existing value
        if routes_for_stop:
            destinations.at[index, 'routes'] = ', '.join(routes_for_stop)

    return destinations


def get_routes_dict(bus_stop):
    # Initialize a dictionary for routes with bus stop details
    routes_dict = defaultdict(list)

    # Iterate through the GeoDataFrame to populate dictionaries
    for index, row in bus_stop.iterrows():
        if pd.isna(row['routes']) or row['routes'] == 'NaN':
            # Add details to the no_route_dict
            routes_dict['no_route'].append(row['sap_name'])
        else:
            # Split routes and add details to the respective route in routes_dict
            routes = str(row['routes']).split(', ')
            for route in routes:
                routes_dict[route].append(
                    row['sap_name'])

    # Convert defaultdict to regular dict for easier use
    routes_dict = dict(routes_dict)
    return routes_dict


def visual_route_allocation(gdf, route_name, mapbox_token):
    """
    Plots a scatter plot on a Mapbox backdrop using Plotly.
    """
    route_dict = get_routes_dict(gdf)
    gdf1 = gdf[gdf['sap_name'].isin(route_dict[route_name])]
    gdf1['route'] = route_name
    gdf2 = gdf[gdf['sap_name'].isin(route_dict['no_route'])]
    gdf2['route'] = 'no route'
    gdf = pd.concat([gdf1, gdf2], ignore_index=True).reset_index(drop=True)
    # Ensure the Mapbox access token is set
    # mapbox_token = get_mapbox_token()
    px.set_mapbox_access_token(mapbox_token)

    # Create the scatter plot
    fig = px.scatter_mapbox(gdf, lat='lat',
                            lon='lon',
                            zoom=10,
                            height=600,
                            width=1000,
                            mapbox_style="light",
                            color='route',
                            hover_data='sap_name'
                            )

    # Show the plot
    fig.show()


def get_route_geometry(G, bus_stop, route_name, city_grid, plot=True):
    routes_dict = get_routes_dict(bus_stop)

    # Filter DataFrame based on specified route name in the 'routes' column
    filtered_df = city_grid[city_grid['routes'].str.contains(route_name)]

    # Extract unique 'sap_names' from the filtered DataFrame
    duplicates_mask = filtered_df.duplicated(subset=['col1'], keep='first')
    bus_stop_nodes = filtered_df[~duplicates_mask]

    # Initialize lists to store the latitude and longitude of the route points
    route_points = []

    # Compute the shortest path between each pair of consecutive bus stops and collect coordinates
    for start_node, end_node in zip(bus_stop_nodes.node[:-1], bus_stop_nodes.node[1:]):
        shortest_path = nx.shortest_path(G, start_node, end_node, weight='length')
        for node in shortest_path:
            point = G.nodes[node]
            route_points.append((point['x'], point['y']))

    # Create a LineString from the route points
    route_geometry = LineString(route_points)

    if plot:
        # Plotting
        fig, ax = plt.subplots(figsize=(5, 5))

        # Plot city grid
        city_grid.plot(ax=ax, color='lightgray', edgecolor='black', alpha=0.5)

        # Plot the route
        gpd.GeoSeries([route_geometry]).plot(ax=ax, linewidth=2, color='blue')

        # Adjust plot limits
        minx, miny, maxx, maxy = city_grid.total_bounds
        ax.set_xlim(minx, maxx)
        ax.set_ylim(miny, maxy)
        ax.set_aspect('equal', adjustable='box')

        # Optionally, customize the plot further (e.g., title)
        ax.set_title(f'Route: {route_name}')

        plt.axis('off')
        plt.show()

    # Return the route geometry for later use
    return route_geometry


def plot_bivariate_choropleth(gdf, route=None, plot_stops=None, color_legends=True):
    gdf = gdf.to_crs(epsg=3857)

    # Plotting the choropleth map
    fig, ax = plt.subplots(figsize=(10, 10))

    # Define specific colors for each category
    colors = {'3A': '#D9CFC6', '2A': '#B7999C', '1A': '#916A71', '3B': '#DB9E9C', '2B': '#B3756D', '1B': '#925C59'}
    alpha = 0.8

    # Plot each category separately
    for category, color in colors.items():
        subset = gdf[gdf['bivariate_cat'] == category]
        subset.plot(ax=ax, color=color, edgecolor='grey', alpha=alpha)

    ax.set_title('Accessibility Index vs Population Density')
    ax.set_axis_off()

    # Plot route if provided
    route_legend_items = []
    if route is not None:
        set1_colors = plt.get_cmap('Set1').colors  # Get colors from 'Set1' colormap
        for idx, (route_name, route_geo) in enumerate(route.items()):
            color = set1_colors[idx % len(set1_colors)]  # Cycle through 'Set1' colors
            route_geo = route_geo.to_crs(epsg=3857)
            route_geo.plot(ax=ax, color=color, label=route_name)
            route_legend_items.append(Line2D([0], [0], color=color, lw=2, label=route_name))

    # Plot busy bus stops if provided
    stop_legend_items = []
    if plot_stops is not None:
        dark2_colors = plt.get_cmap('Dark2').colors  # Get colors from 'Dark2' colormap
        for idx, (stops_label, stops_geo) in enumerate(plot_stops.items()):
            color = dark2_colors[idx % len(dark2_colors)]  # Cycle through 'Dark2' colors
            stops_geo = stops_geo.to_crs(epsg=3857)
            stops_geo.plot(ax=ax, color=color, markersize=50, label=stops_label)
            stop_legend_items.append(Line2D([0], [0], color=color, marker='o', linestyle='None', markersize=10, label=stops_label))

            # Annotate bus stop names
            for x, y, sap_name in zip(stops_geo.geometry.x, stops_geo.geometry.y, stops_geo['sap_name']):
                ax.text(x, y, sap_name, fontsize=8, ha='right', va='bottom')

    # Combine all legend items
    all_legend_items = route_legend_items + stop_legend_items

    # Create the legend with custom items
    ax.legend(handles=all_legend_items, loc='lower right')

    # Add a basemap
    provider = ctx.providers['OpenStreetMap']['Mapnik']
    ctx.add_basemap(ax, source=provider)

    if color_legends:
        img2 = fig  # refer to the main figure
        ax2 = fig.add_axes([0.24, 0.15, 0.1, 0.1])  # add new axes to place the legend there

        # Column 1
        ax2.axvspan(xmin=0, xmax=0.33, ymin=0, ymax=0.33, alpha=alpha, color=colors['1A'])
        ax2.axvspan(xmin=0, xmax=0.33, ymin=0.33, ymax=0.66, alpha=alpha, color=colors['2A'])
        ax2.axvspan(xmin=0, xmax=0.33, ymin=0.66, ymax=1, alpha=alpha, color=colors['3A'])

        # Column 2
        ax2.axvspan(xmin=0.33, xmax=0.66, ymin=0, ymax=0.33, alpha=alpha, color=colors['1B'])
        ax2.axvspan(xmin=0.33, xmax=0.66, ymin=0.33, ymax=0.66, alpha=alpha, color=colors['2B'])
        ax2.axvspan(xmin=0.33, xmax=0.66, ymin=0.66, ymax=1, alpha=alpha, color=colors['3B'])

        # Step 3: annotate the legend
        ax2.tick_params(axis='both', which='both', length=0, )  # remove ticks from the big box
        ax2.annotate("", xy=(0, 1), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="black"))  # draw arrow for x
        ax2.annotate("", xy=(1, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="black"))  # draw arrow for y
        ax2.text(s='Population/Sq.km', x=0.1, y=-0.25)  # annotate x-axis
        ax2.text(s='PTAL', x=-0.25, y=0.1, rotation=90)  # annotate y-axis
        ax2.set_axis_off()
    plt.show()
