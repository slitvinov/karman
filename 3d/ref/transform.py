import math
import re
import sys

# Re = 2120
t = 100 / 47, -10 / 47, 0
r = 0, 0, -90
s = 2 / 47, 2 / 47, 2 / 47
p = math.radians(r[2])
T = [[math.cos(p), -math.sin(p)], [math.sin(p), math.cos(p)]]

for path in sys.argv[1:]:
    with open(path, "r") as f:
        title = re.sub('TITLE[ \t]*=[\ t]*"', '', f.readline())
        title = re.sub('"[ \t]*\n', '', title)
        variables = re.sub('VARIABLES[ \t]=[ \t]*', '', f.readline())
        variables = re.sub('[ \t]*\n', '', variables)
        variables = re.split(",[ \t]*", variables)
        variables = [re.sub('(^")|("$)', "", x) for x in variables]
        zone = re.sub('ZONE[ \t]*', '', f.readline())
        zone = re.sub('[ \t]*\n', '', zone)
        zone = re.split(",[ \t]*", zone)
        zone = dict(x.split("=") for x in zone)
        nx = int(zone["I"])
        ny = int(zone["J"])
        nz = int(zone["K"])
        data = [[] for x in variables]
        for line in f:
            for d, x in zip(data, line.split()):
                x = float(x)
                d.append(x)
        assert len(data[0]) == nx * ny
        Data = dict(zip(variables, data))

    # rotate
    Data["x"], Data["y"] = zip(*[(T[0][0] * x + T[0][1] * y,
                                  T[1][0] * x + T[1][1] * y)
                                 for x, y in zip(Data["x"], Data["y"])])
    # scale and translate
    Data["x"] = [s[0] * r + t[0] for r in Data["x"]]
    Data["y"] = [s[1] * r + t[1] for r in Data["y"]]
    Data["z"] = [s[2] * r + t[2] for r in Data["z"]]

    Data["Vx"], Data["Vy"] = Data["Vy"], Data["Vx"]

    xdmf_path = re.sub("\.dat$", "", path) + ".xdmf2"
    with open(xdmf_path, "w") as f:
        f.write("""\
    <Xdmf
        Version="2">
      <Domain>
        <Grid>
          <Topology
              TopologyType="3DSMesh"
              Dimensions="%ld %ld %ld"/>
          <Geometry
              GeometryType="XYZ">
            <DataItem
                Dimensions="%ld 3">
    """ % (nz, ny, nx, nx * ny * nz))
        for x, y, z in zip(Data["x"], Data["y"], Data["z"]):
            f.write("""\
              %.16e %.16e %.16e
    """ % (x, y, z))
        f.write("""\
            </DataItem>
          </Geometry>
    """)
        f.write("""\
          <Attribute
              Name="u"
              AttributeType="Vector">
            <DataItem
                Dimensions="1 %ld %ld 3">
    """ % (ny, nx))
        for x, y, z in zip(Data["Vx"], Data["Vy"], Data["Vz"]):
            f.write("""\
            %.16e %.16e %.16e
    """ % (x, y, z))
        f.write("""\
            </DataItem>
          </Attribute>
        </Grid>
      </Domain>
    </Xdmf>
    """)
