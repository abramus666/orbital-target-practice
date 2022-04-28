
import math, random, time
import tkinter as tk
import tkinter.ttk as ttk

FRAME_DELAY_MS      = 10

ROCK_COLOR          = (0.5, 0.4, 0.3)
HIGHLIGHT_COLOR     = (0.25, 0.25, 1.0)

MESSAGE_DURATION    = 3.0
TARGET_AREA_SIZE    = 4.0
DAMAGE_MIN          = 5
DAMAGE_MAX          = 50

SPAWN_RATE_DURATION = 60.0
SPAWN_RATES         = (0.8, 0.7, 0.6, 0.5, 0.4, 0.3)

VIEW_SHAKE_SCALE    = 0.1
VIEW_SHAKE_DURATION = 1.0

ROCK_SPAWN_DISTANCE = 200.0
ROCK_RADIUS         = 1.0
ROCK_SPEED          = 20.0
ROCK_ANGULAR_SPEED  = 2.0

BOOM_FRAG_COUNT     = 10
BOOM_STEP_COUNT     = 10
BOOM_RADIUS         = 0.5
BOOM_EXPANSION_TIME = 0.2
BOOM_DURATION       = 0.8

BEAM_DURATION       = 0.1
BEAM_START_POSITION = (0.0, -0.6, 1.0)
BEAM_WIDTH          = 0.1
BEAM_STRIP_COUNT    = 15
BEAM_MAX_COLOR      = (0.5, 1.0, 1.0)
BEAM_MIN_COLOR      = (0.0, 0.5, 1.0)

LOD0_THRESHOLD_PX   = 5
LOD1_THRESHOLD_PX   = 30
LOD2_THRESHOLD_PX   = 90
LOD_COUNT           = 3

DEFORM_COUNT        = 10
DEFORM_MAJOR_CHOICE = (0.8, 1.25)
DEFORM_EXTRA_RANGE  = (0.98, 1.02)

#===============================================================================
class Vertex:

   def __init__(self, x, y, z):
      self.sphere_coordinates = (x,y,z)

   def applyDeformations(self, major_deforms, extra_deform):
      x,y,z = self.sphere_coordinates
      cumulative_deform = 1.0
      for deform in major_deforms:
         dx = x - deform[0]
         dy = y - deform[1]
         dz = z - deform[2]
         dist_sqr = dx*dx + dy*dy + dz*dz
         strength = 1.0 / (1.0 + 5.0 * dist_sqr * dist_sqr)
         cumulative_deform *= strength * (deform[3] - 1.0) + 1.0
      self.x = x * cumulative_deform * extra_deform
      self.y = y * cumulative_deform * extra_deform
      self.z = z * cumulative_deform * extra_deform

#===============================================================================
class Triangle:

   def __init__(self, v1, v2, v3):
      self.v1 = v1
      self.v2 = v2
      self.v3 = v3

#===============================================================================
class Sphere:

   def __init__(self, radius):
      self.vertices = []
      self.vertex_counts = []
      self.triangles = []
      vertex_indexes = [[],[],[],[],[],[],[]]
      average_indexes = {}

      # Rotate around the Y axis.
      def rotateY(vertex, angle):
         c = math.cos(angle)
         s = math.sin(angle)
         x = c*vertex[0] + s*vertex[2]
         y = vertex[1]
         z = -s*vertex[0] + c*vertex[2]
         return (x,y,z)

      # Rotate around the Z axis.
      def rotateZ(vertex, angle):
         c = math.cos(angle)
         s = math.sin(angle)
         x = c*vertex[0] - s*vertex[1]
         y = s*vertex[0] + c*vertex[1]
         z = vertex[2]
         return (x,y,z)

      # Add vertex to the array, and return
      # index of the added vertex within the array.
      def addVertex(vertex):
         ix = len(self.vertices)
         self.vertices.append(vertex)
         return ix

      # Calculate coordinates of a vertex located
      # on a sphere between two given vertices.
      def makeAverageVertex(v1, v2):
         x = (v1[0] + v2[0]) / 2.0
         y = (v1[1] + v2[1]) / 2.0
         z = (v1[2] + v2[2]) / 2.0
         length = math.sqrt(x*x + y*y + z*z)
         x *= (radius / length)
         y *= (radius / length)
         z *= (radius / length)
         return (x,y,z)

      # Add vertex located on a sphere between two vertices to the
      # array, unless this vertex is already present in the array.
      # Return index of this vertex within the array.
      def addAverageVertex(v1_ix, v2_ix):
         if v1_ix > v2_ix:
            v1_ix, v2_ix = v2_ix, v1_ix
         key = '{}.{}'.format(v1_ix, v2_ix)
         if key not in average_indexes:
            average_indexes[key] = addVertex(makeAverageVertex(self.vertices[v1_ix], self.vertices[v2_ix]))
         return average_indexes[key]

      # For level 0, simply add a triangle defined by the given vertices to the array.
      # For levels above 0, add 4 triangles for each triangle added on the previous
      # level to achieve higher level of detail, and add missing vertices.
      def addTriangles(level, v1_ix, v2_ix, v3_ix):
         if level > 0:
            v12_ix = addAverageVertex(v1_ix, v2_ix)
            v23_ix = addAverageVertex(v2_ix, v3_ix)
            v31_ix = addAverageVertex(v3_ix, v1_ix)
            addTriangles(level-1, v1_ix,  v12_ix, v31_ix)
            addTriangles(level-1, v12_ix, v2_ix,  v23_ix)
            addTriangles(level-1, v31_ix, v23_ix, v3_ix)
            addTriangles(level-1, v12_ix, v23_ix, v31_ix)
         else:
            self.triangles[-1].append((v1_ix, v2_ix, v3_ix))

      # Add top and bottom vertices.
      vertex_indexes[0].append(addVertex((0, radius,0)))
      vertex_indexes[6].append(addVertex((0,-radius,0)))
      # Add remaining vertices for level 0 sphere.
      for vix,cnt,degz,degy in ((1,5,60,0), (2,10,30,0), (3,10,0,-18), (4,10,-30,0), (5,5,-60,0)):
         for i in range(cnt):
            v = rotateY(rotateZ((radius,0,0), degz*math.pi/180), degy*math.pi/180 + i*2*math.pi/cnt)
            vertex_indexes[vix].append(addVertex(v))
      # Add triangles for all levels of detail, starting with 0.
      # For levels above 0, missing vertices are added when necessary.
      # Note that vertices for level 0 are shared between triangles for all levels,
      # vertices for level 1 are shared between triangles for level 1 and above,
      # and so on.
      for level in range(LOD_COUNT):
         # Add empty array for triangles for the current level.
         self.triangles.append([])
         # Add triangles for the current level, adding missing vertices when necessary.
         for i in range(5):
            # 1st row (top).
            addTriangles(level, vertex_indexes[0][0],       vertex_indexes[1][i],          vertex_indexes[1][(i+1)%5])
            # 2nd row.
            addTriangles(level, vertex_indexes[1][i],       vertex_indexes[2][i*2],        vertex_indexes[2][i*2+1])
            addTriangles(level, vertex_indexes[1][i],       vertex_indexes[2][i*2+1],      vertex_indexes[1][(i+1)%5])
            addTriangles(level, vertex_indexes[1][(i+1)%5], vertex_indexes[2][i*2+1],      vertex_indexes[2][(i*2+2)%10])
            # 3rd row.
            addTriangles(level, vertex_indexes[2][i*2],     vertex_indexes[3][i*2],        vertex_indexes[3][i*2+1])
            addTriangles(level, vertex_indexes[2][i*2],     vertex_indexes[3][i*2+1],      vertex_indexes[2][i*2+1])
            addTriangles(level, vertex_indexes[2][i*2+1],   vertex_indexes[3][i*2+1],      vertex_indexes[3][(i*2+2)%10])
            addTriangles(level, vertex_indexes[2][i*2+1],   vertex_indexes[3][(i*2+2)%10], vertex_indexes[2][(i*2+2)%10])
            # 4th row.
            addTriangles(level, vertex_indexes[4][i*2],     vertex_indexes[3][i*2+1],      vertex_indexes[3][i*2])
            addTriangles(level, vertex_indexes[4][i*2],     vertex_indexes[4][i*2+1],      vertex_indexes[3][i*2+1])
            addTriangles(level, vertex_indexes[4][i*2+1],   vertex_indexes[3][(i*2+2)%10], vertex_indexes[3][i*2+1])
            addTriangles(level, vertex_indexes[4][i*2+1],   vertex_indexes[4][(i*2+2)%10], vertex_indexes[3][(i*2+2)%10])
            # 5th row.
            addTriangles(level, vertex_indexes[5][i],       vertex_indexes[4][i*2+1],      vertex_indexes[4][i*2])
            addTriangles(level, vertex_indexes[5][i],       vertex_indexes[5][(i+1)%5],    vertex_indexes[4][i*2+1])
            addTriangles(level, vertex_indexes[5][(i+1)%5], vertex_indexes[4][(i*2+2)%10], vertex_indexes[4][i*2+1])
            # 6th row (bottom).
            addTriangles(level, vertex_indexes[6][0],       vertex_indexes[5][(i+1)%5],    vertex_indexes[5][i])
         # Save number of vertices for each level of 
         self.vertex_counts.append(len(self.vertices))

#===============================================================================
class Rock:

   def __init__(self, sphere):
      # Vertices are shared between triangles for all levels of detail.
      self.vertices = [Vertex(*vert) for i,vert in enumerate(sphere.vertices)]
      self.vertex_counts = sphere.vertex_counts
      # For each level of detail, create array of triangles.
      self.triangles = []
      for triangles in sphere.triangles:
         self.triangles.append([Triangle(self.vertices[tri[0]], self.vertices[tri[1]], self.vertices[tri[2]]) for tri in triangles])
      # Initialize remaining data.
      self.major_deforms = [[0,0,0,0] for i in range(DEFORM_COUNT)]
      self.alive = False
      self.angle = [0,0,0]
      self.rotation = [0,0,0]
      self.position = [0,0,0]
      self.target = [0,0,0]
      self.velocity = [0,0,0]

   def destroy(self):
      self.alive = False

   def select(self):
      self.color = HIGHLIGHT_COLOR

   def deselect(self):
      self.color = ROCK_COLOR

   def spawn(self):
      self.alive = True
      self.color = ROCK_COLOR
      # The base shape of a rock is a sphere, that is deformed.
      # Major deformations are to change the overall shape of a sphere.
      for deform in self.major_deforms:
         ix = random.randrange(len(self.vertices))
         deform[0] = self.vertices[ix].sphere_coordinates[0]
         deform[1] = self.vertices[ix].sphere_coordinates[1]
         deform[2] = self.vertices[ix].sphere_coordinates[2]
         deform[3] = random.choice(DEFORM_MAJOR_CHOICE)
      # Extra deformations are to add noise and make a rock less regular.
      for vert in self.vertices:
         extra_deform = random.uniform(*DEFORM_EXTRA_RANGE)
         vert.applyDeformations(self.major_deforms, extra_deform)
      # Randomly choose rotation and initial position.
      self.rotation[0] = random.uniform(-ROCK_ANGULAR_SPEED, ROCK_ANGULAR_SPEED)
      self.rotation[1] = random.uniform(-ROCK_ANGULAR_SPEED, ROCK_ANGULAR_SPEED)
      self.rotation[2] = random.uniform(-ROCK_ANGULAR_SPEED, ROCK_ANGULAR_SPEED)
      self.position[0] = random.uniform(-ROCK_SPAWN_DISTANCE*0.4, ROCK_SPAWN_DISTANCE*0.4)
      self.position[1] = random.uniform(-ROCK_SPAWN_DISTANCE*0.4, ROCK_SPAWN_DISTANCE*0.4)
      self.position[2] = ROCK_SPAWN_DISTANCE
      # Randomly choose target point, and calculate velocity vector for it.
      self.target[0] = random.uniform(-TARGET_AREA_SIZE/2, TARGET_AREA_SIZE/2)
      self.target[1] = random.uniform(-TARGET_AREA_SIZE/2, TARGET_AREA_SIZE/2)
      self.target[2] = 0
      dx = self.target[0] - self.position[0]
      dy = self.target[1] - self.position[1]
      dz = self.target[2] - self.position[2]
      length = math.sqrt(dx*dx + dy*dy + dz*dz)
      self.velocity[0] = ROCK_SPEED * dx/length
      self.velocity[1] = ROCK_SPEED * dy/length
      self.velocity[2] = ROCK_SPEED * dz/length

   def update(self, dt):
      self.angle[0] += self.rotation[0] * dt
      self.angle[1] += self.rotation[1] * dt
      self.angle[2] += self.rotation[2] * dt
      self.position[0] += self.velocity[0] * dt
      self.position[1] += self.velocity[1] * dt
      self.position[2] += self.velocity[2] * dt
      if self.position[2] <= 0:
         # Crashed into the ground.
         self.alive = False

   def renderInPanel_simple(self, panel):
      # Weak perspective projection.
      x = self.position[0] / self.position[2]
      y = self.position[1] / self.position[2]
      r = ROCK_RADIUS / self.position[2]
      # Convert to view coordinates.
      view_x = int((0.5 + x) * panel.view_scale + panel.view_x)
      view_y = int((0.5 - y) * panel.view_scale + panel.view_y)
      view_r = int(r * panel.view_scale)
      # Check if this rock is near enough to be seen.
      if view_r > 0:
         # Determine color.
         brightness = 0.75
         r = int(self.color[0] * brightness * 255.0)
         g = int(self.color[1] * brightness * 255.0)
         b = int(self.color[2] * brightness * 255.0)
         hex_color = '#{:02x}{:02x}{:02x}'.format(r,g,b)
         # Get oval.
         oval_ix = panel.allocateOvals(1)
         oval_id = panel.oval_pool[oval_ix]
         # Update oval coordinates and color.
         panel.canvas.coords(oval_id, view_x-view_r, view_y-view_r, view_x+view_r, view_y+view_r)
         panel.canvas.itemconfig(oval_id, state = 'normal', fill = hex_color)
         panel.canvas.tag_lower(oval_id)

   def renderInPanel(self, panel):
      # Level of detail depends on the distance and view panel size.
      view_radius = panel.view_scale / self.position[2]
      if view_radius < LOD0_THRESHOLD_PX:
         self.renderInPanel_simple(panel)
         return
      elif view_radius < LOD1_THRESHOLD_PX:
         array_ix = 0
      elif view_radius < LOD2_THRESHOLD_PX:
         array_ix = 1
      else:
         array_ix = LOD_COUNT-1
      # Precompute cos/sin values.
      cx = math.cos(self.angle[0])
      cy = math.cos(self.angle[1])
      cz = math.cos(self.angle[2])
      sx = math.sin(self.angle[0])
      sy = math.sin(self.angle[1])
      sz = math.sin(self.angle[2])
      # Process vertices.
      for vert in self.vertices[:self.vertex_counts[array_ix]]:
         # Rotate around the X axis.
         y_ = cx*vert.y - sx*vert.z
         z_ = sx*vert.y + cx*vert.z
         # Rotate around the Y axis.
         x_ = cy*vert.x + sy*z_
         z = -sy*vert.x + cy*z_
         # Rotate around the Z axis.
         x = cz*x_ - sz*y_
         y = sz*x_ + cz*y_
         # Move.
         vert.camera_x = x + self.position[0]
         vert.camera_y = y + self.position[1]
         vert.camera_z = z + self.position[2]
         # Weak perspective projection.
         x = vert.camera_x / vert.camera_z
         y = vert.camera_y / vert.camera_z
         # Convert to view coordinates.
         vert.view_x = int((0.5 + x) * panel.view_scale + panel.view_x)
         vert.view_y = int((0.5 - y) * panel.view_scale + panel.view_y)
         vert.visible = (
            (vert.view_x > 0) and (vert.view_x < panel.width) and
            (vert.view_y > 0) and (vert.view_y < panel.height))
      # Process triangles.
      num_visible_triangles = 0
      for tri in self.triangles[array_ix]:
         tri.distance_z = (tri.v1.camera_z + tri.v2.camera_z + tri.v3.camera_z) / 3.0
         tri.visible = (tri.distance_z > 0) and (tri.v1.visible or tri.v2.visible or tri.v3.visible)
         if tri.visible:
            # 1st triangle vector.
            ax = tri.v2.camera_x - tri.v1.camera_x
            ay = tri.v2.camera_y - tri.v1.camera_y
            az = tri.v2.camera_z - tri.v1.camera_z
            # 2nd triangle vector.
            bx = tri.v3.camera_x - tri.v1.camera_x
            by = tri.v3.camera_y - tri.v1.camera_y
            bz = tri.v3.camera_z - tri.v1.camera_z
            # Cross product to determine normal.
            x = ay*bz - az*by
            y = az*bx - ax*bz
            z = ax*by - ay*bx
            # Since camera always looks along the Z axis, we don't need dot product here.
            # Simply use Z coordinate of the unit normal to determine brightness.
            length = math.sqrt(x*x + y*y + z*z)
            tri.brightness = -z / length
            if tri.brightness > 0:
               num_visible_triangles += 1
            else:
               tri.visible = False
      # Render triangles.
      poly_ix = panel.allocatePolygons(num_visible_triangles)
      for tri in self.triangles[array_ix]:
         if tri.visible:
            # Determine color.
            r = int(self.color[0] * tri.brightness * 255.0)
            g = int(self.color[1] * tri.brightness * 255.0)
            b = int(self.color[2] * tri.brightness * 255.0)
            hex_color = '#{:02x}{:02x}{:02x}'.format(r,g,b)
            # Get polygon.
            tri.poly_id = panel.polygon_pool[poly_ix]
            poly_ix += 1
            # Update polygon coordinates and color.
            panel.canvas.coords(tri.poly_id, tri.v1.view_x, tri.v1.view_y, tri.v2.view_x, tri.v2.view_y, tri.v3.view_x, tri.v3.view_y)
            panel.canvas.itemconfig(tri.poly_id, state = 'normal', fill = hex_color)
      # Sort triangles by their distance, so that they are rendered in the correct order.
      self.triangles[array_ix].sort(key = lambda tri: tri.distance_z)
      for tri in self.triangles[array_ix]:
         if tri.visible:
            panel.canvas.tag_lower(tri.poly_id)

#===============================================================================
class Boom:

   def __init__(self):
      self.alive = False
      self.age = 0
      self.position = [0,0,0]
      self.velocity = [0,0,0]
      self.sub_elements = [[0,0,0,0,0,0,0,0] for i in range(BOOM_FRAG_COUNT * BOOM_STEP_COUNT)]

   def spawn(self, position, velocity):
      self.alive = True
      self.age = 0
      self.position[0] = position[0]
      self.position[1] = position[1]
      self.position[2] = position[2]
      self.velocity[0] = velocity[0] / 10.0
      self.velocity[1] = velocity[1] / 10.0
      self.velocity[2] = velocity[2] / 10.0
      # Delta time.
      dt = BOOM_EXPANSION_TIME / BOOM_STEP_COUNT
      # Create fragments.
      for frag_ix in range(BOOM_FRAG_COUNT):
         x,y,z = position
         # Randomly choose direction for the successive steps of this fragment.
         vx = velocity[0] + random.uniform(-ROCK_SPEED, ROCK_SPEED)
         vy = velocity[1] + random.uniform(-ROCK_SPEED, ROCK_SPEED)
         vz = velocity[2]
         length = math.sqrt(vx*vx + vy*vy + vz*vz)
         vx *= ROCK_SPEED / length
         vy *= ROCK_SPEED / length
         vz *= ROCK_SPEED / length
         # Determine properties for the successive steps of this fragment.
         for step_ix in range(BOOM_STEP_COUNT):
            # X/Y changes fade away for later steps.
            x += vx * dt * (1.0 - step_ix/BOOM_STEP_COUNT)
            y += vy * dt * (1.0 - step_ix/BOOM_STEP_COUNT)
            z += vz * dt
            sub = self.sub_elements[BOOM_STEP_COUNT * frag_ix + step_ix]
            # Position of the current subelement.
            sub[0] = x
            sub[1] = y
            sub[2] = z
            # Velocity of the current subelement. Add some randomness in X/Y plane.
            sub[3] = self.velocity[0] + random.uniform(-BOOM_RADIUS/2, BOOM_RADIUS/2)
            sub[4] = self.velocity[1] + random.uniform(-BOOM_RADIUS/2, BOOM_RADIUS/2)
            sub[5] = self.velocity[2]
            # Delay of the current subelement. Range [0, BOOM_EXPANSION_TIME].
            sub[6] = BOOM_EXPANSION_TIME * step_ix / (BOOM_STEP_COUNT - 1.0)
            # Radius of the current subelement. Range [1, 0).
            sub[7] = BOOM_RADIUS * math.sqrt((BOOM_STEP_COUNT - step_ix) / BOOM_STEP_COUNT)
      # Sort subelements by their distance, so that they are rendered in the correct order.
      self.sub_elements.sort(key = lambda sub: sub[2])

   def update(self, dt):
      # Update position (including all subelements).
      self.position[0] += self.velocity[0] * dt
      self.position[1] += self.velocity[1] * dt
      self.position[2] += self.velocity[2] * dt
      for sub in self.sub_elements:
         sub[0] += sub[3] * dt
         sub[1] += sub[4] * dt
         sub[2] += sub[5] * dt
      # Update age and check whether we are still alive.
      self.age += dt
      if self.age > (BOOM_EXPANSION_TIME + BOOM_DURATION):
         self.alive = False

   def renderInPanel(self, panel):
      for sub in self.sub_elements:
         delay = sub[6]
         radius = sub[7]
         relative_age = (self.age - delay) / BOOM_DURATION
         # Skip this subelement if it should not be active now.
         if (relative_age >= 0.0) and (relative_age <= 1.0) and (sub[2] > 0.0):
            # Determine current size.
            size = relative_age
            for i in range(4):
               size = math.sin(size * math.pi/2)
            radius *= size
            # Weak perspective projection.
            x = sub[0] / sub[2]
            y = sub[1] / sub[2]
            radius /= sub[2]
            # Convert to view coordinates.
            view_x = int((0.5 + x) * panel.view_scale + panel.view_x)
            view_y = int((0.5 - y) * panel.view_scale + panel.view_y)
            view_r = int(radius * panel.view_scale)
            # Check if this subelement is near enough to be seen.
            if view_r > 0:
               # Determine current color.
               r = int(min(max(3.0 - 3.0 * relative_age, 0.0), 1.0) * 255.0)
               g = int(min(max(2.0 - 3.0 * relative_age, 0.0), 1.0) * 255.0)
               b = int(min(max(1.5 - 3.0 * relative_age, 0.0), 1.0) * 255.0)
               hex_color = '#{:02x}{:02x}{:02x}'.format(r,g,b)
               # Get oval.
               oval_ix = panel.allocateOvals(1)
               oval_id = panel.oval_pool[oval_ix]
               # Update oval coordinates and color.
               panel.canvas.coords(oval_id, view_x-view_r, view_y-view_r, view_x+view_r, view_y+view_r)
               panel.canvas.itemconfig(oval_id, state = 'normal', fill = hex_color)
               panel.canvas.tag_lower(oval_id)

#===============================================================================
class Beam:

   def __init__(self):
      self.alive = False
      self.age = 0
      self.position = [0,0,0]
      self.coordinates = [[0,0] for i in range(BEAM_STRIP_COUNT+2)]

   def spawn(self, target):
      self.alive = True
      self.age = 0
      self.position[0] = target[0]
      self.position[1] = target[1]
      self.position[2] = target[2]

   def update(self, dt):
      self.age += dt
      if self.age > BEAM_DURATION:
         self.alive = False

   def renderInPanel(self, panel):
      # Unit vector directed at the target.
      dx = self.position[0] - BEAM_START_POSITION[0]
      dy = self.position[1] - BEAM_START_POSITION[1]
      dz = self.position[2] - BEAM_START_POSITION[2]
      length = math.sqrt(dx*dx + dy*dy + dz*dz)
      dx /= length
      dy /= length
      dz /= length
      # Weak perspective projection.
      # Target position is last.
      self.coordinates[-1][0] = self.position[0] / self.position[2]
      self.coordinates[-1][1] = self.position[1] / self.position[2]
      # Multiple start positions to render strips of different colors.
      for coord_ix in range(BEAM_STRIP_COUNT+1):
         shift = (float(coord_ix) / BEAM_STRIP_COUNT) * 2.0 - 1.0 # Range [-1,1].
         shift = math.copysign(1, shift) * math.sqrt(abs(shift)) # Same range but not linear.
         # Calculate points by shifting along the vector that is perpendicular
         # to both the vector directed at target and the Y axis.
         x = BEAM_START_POSITION[0] - dz * shift * BEAM_WIDTH/2
         y = BEAM_START_POSITION[1]
         z = BEAM_START_POSITION[2] + dx * shift * BEAM_WIDTH/2
         self.coordinates[coord_ix][0] = x / z
         self.coordinates[coord_ix][1] = y / z
      # Convert to view coordinates.
      for coordinate in self.coordinates:
         x,y = coordinate
         coordinate[0] = int((0.5 + x) * panel.view_scale + panel.view_x)
         coordinate[1] = int((0.5 - y) * panel.view_scale + panel.view_y)
      # Render strips.
      poly_ix = panel.allocatePolygons(BEAM_STRIP_COUNT)
      for coord_ix in range(BEAM_STRIP_COUNT):
         # Determine color.
         factor = (float(coord_ix) / (BEAM_STRIP_COUNT-1)) * 2.0 - 1.0 # Range [-1,1].
         factor = abs(factor) # 0 for middle index, 1 for boundary indexes.
         r = int((BEAM_MAX_COLOR[0] * (1.0 - factor) + BEAM_MIN_COLOR[0] * factor) * 255.0)
         g = int((BEAM_MAX_COLOR[1] * (1.0 - factor) + BEAM_MIN_COLOR[1] * factor) * 255.0)
         b = int((BEAM_MAX_COLOR[2] * (1.0 - factor) + BEAM_MIN_COLOR[2] * factor) * 255.0)
         hex_color = '#{:02x}{:02x}{:02x}'.format(r,g,b)
         # Get polygon.
         poly_id = panel.polygon_pool[poly_ix]
         poly_ix += 1
         # Update polygon coordinates and color.
         x0,y0 = self.coordinates[-1]
         x1,y1 = self.coordinates[coord_ix]
         x2,y2 = self.coordinates[coord_ix+1]
         panel.canvas.coords(poly_id, x0, y0, x1, y1, x2, y2)
         panel.canvas.itemconfig(poly_id, state = 'normal', fill = hex_color)
         panel.canvas.tag_lower(poly_id)

#===============================================================================
class ViewPanel:

   def __init__(self, parent, row, column, width, height):
      # Event related.
      self.mouse_position = [0,0]
      self.mouse_pressed = False
      self.start_pressed = False
      # Time related.
      self.shake_scale = 0
      self.text_delay  = 0
      # Allocation.
      self.polygon_pool = []
      self.oval_pool = []
      self.num_polygons = 0
      self.num_ovals = 0
      # Frame.
      self.frame = ttk.Frame(parent)
      self.frame.grid(row = row, column = column, sticky = 'NSEW')
      self.frame.columnconfigure(0, weight = 1)
      self.frame.rowconfigure(0, weight = 1)
      # Canvas.
      self.canvas = tk.Canvas(self.frame, width = width, height = height, cursor = 'plus', background = 'black')
      self.canvas.grid(sticky = 'NSEW')
      self.canvas.bind('<Configure>', self.resizeEvent)
      self.canvas.bind('<Motion>', self.mouseMotionEvent)
      self.canvas.bind('<Button-1>', self.mousePressEvent)
      self.canvas.bind_all('<Key-space>', self.startPressEvent)
      # HUD.
      self.text_top    = self.canvas.create_text(0,0, anchor = tk.N,      state = 'hidden', font = 'consolas 12', fill = 'white')
      self.text_center = self.canvas.create_text(0,0, anchor = tk.CENTER, state = 'hidden', font = 'consolas 18', fill = 'white')
      self.text_sw     = self.canvas.create_text(0,0, anchor = tk.SW,     state = 'hidden', font = 'consolas 12', fill = 'white')
      self.text_se     = self.canvas.create_text(0,0, anchor = tk.SE,     state = 'hidden', font = 'consolas 12', fill = 'white')
      # View.
      self.setupView(width, height)

   #----------------------------------------------------------------------------

   def setupView(self, width, height):
      self.width = width
      self.height = height
      # Prepare parameters for conversion to view coordinates.
      # X and Y in range [-0.5, 0.5] are always in view,
      # others may or may not depending on the view ratio.
      if width > height:
         self.view_scale = height
         self.view_origin_x = int((width - height) / 2.0)
         self.view_origin_y = 0
      else:
         self.view_scale = width
         self.view_origin_x = 0
         self.view_origin_y = int((height - width) / 2.0)
      # Need copies for view shaking.
      self.view_x = self.view_origin_x
      self.view_y = self.view_origin_y
      # Update text positions.
      self.canvas.coords(self.text_top,    width/2,  10)
      self.canvas.coords(self.text_center, width/2,  height/2)
      self.canvas.coords(self.text_sw,     10,       height-10)
      self.canvas.coords(self.text_se,     width-10, height-10)

   def resizeEvent(self, event):
      self.setupView(event.width, event.height)

   def mouseMotionEvent(self, event):
      x = event.x - self.view_x
      y = event.y - self.view_y
      # Convert from view coordinates.
      self.mouse_position[0] = float(x) / float(self.view_scale) - 0.5
      self.mouse_position[1] = 0.5 - float(y) / float(self.view_scale)

   def mousePressEvent(self, event):
      self.mouseMotionEvent(event)
      self.mouse_pressed = True

   def startPressEvent(self, event):
      self.start_pressed = True

   #----------------------------------------------------------------------------
   def update(self, dt):
      if self.shake_scale > 0:
         self.shake_scale = max(self.shake_scale - (VIEW_SHAKE_SCALE / VIEW_SHAKE_DURATION) * dt, 0)
         dx = self.shake_scale * self.view_scale
         dy = self.shake_scale * self.view_scale
         self.view_x = self.view_origin_x + random.uniform(-dx, dx)
         self.view_y = self.view_origin_y + random.uniform(-dy, dy)
      if self.text_delay > 0:
         self.text_delay -= dt
         if not (self.text_delay > 0):
            self.canvas.itemconfig(self.text_top, state = 'hidden')

   def shakeView(self, scale):
      if self.shake_scale < scale:
         self.shake_scale = scale

   def printMessage(self, message):
      self.text_delay = MESSAGE_DURATION
      self.canvas.itemconfig(self.text_top, state = 'normal', text = message)

   def showStartView(self):
      self.canvas.itemconfig(self.text_center, state = 'normal', text = 'Press SPACE to start')

   def hideStartView(self):
      self.canvas.itemconfig(self.text_center, state = 'hidden')

   def setScore(self, score):
      self.canvas.itemconfig(self.text_sw, state = 'normal', text = 'Score: {:,}'.format(score))

   def setHealth(self, health):
      self.canvas.itemconfig(self.text_se, state = 'normal', text = 'Health: {}'.format(health))

   #----------------------------------------------------------------------------
   def allocatePolygons(self, count):
      ix = self.num_polygons
      self.num_polygons += count
      while len(self.polygon_pool) < self.num_polygons:
         self.polygon_pool.append(self.canvas.create_polygon(0,0,0,0,0,0))
      return ix

   def allocateOvals(self, count):
      ix = self.num_ovals
      self.num_ovals += count
      while len(self.oval_pool) < self.num_ovals:
         self.oval_pool.append(self.canvas.create_oval(0,0,0,0, width = 0))
      return ix

   def deallocateAll(self):
      self.num_polygons = 0
      self.num_ovals = 0

   def noMoreAllocations(self):
      for i in range(self.num_polygons, len(self.polygon_pool)):
         self.canvas.itemconfig(self.polygon_pool[i], state = 'hidden')
      for i in range(self.num_ovals, len(self.oval_pool)):
         self.canvas.itemconfig(self.oval_pool[i], state = 'hidden')

#===============================================================================
class GameLogic:

   def __init__(self, panel):
      self.panel = panel
      self.sphere = Sphere(ROCK_RADIUS)
      self.rocks = []
      self.booms = []
      self.beam = Beam()
      self.rendered_objects = [self.beam]
      self.current_time = time.perf_counter()
      self.resetGame()

   def resetGame(self):
      # Kill all objects.
      for rock in self.rocks:
         rock.alive = False
      for boom in self.booms:
         boom.alive = False
      self.beam.alive = False
      # Reset variables.
      self.start_time = time.perf_counter()
      self.spawn_delay = 0
      self.consecutive_hits = 0
      self.hit_queue = [None,None,None,None]
      self.score_queue = []
      self.score_per_kill = 10
      self.score = 0
      self.health = 100
      self.panel.setScore(self.score)
      self.panel.setHealth(self.health)

   def runFrame(self):
      # Reset variables.
      self.selected_rock = None
      self.unused_rock = None
      self.unused_boom = None
      # Handle time.
      self.previous_time = self.current_time
      self.current_time = time.perf_counter()
      dt = self.current_time - self.previous_time
      # Game restart.
      if self.panel.start_pressed and not (self.health > 0):
         self.panel.hideStartView()
         self.resetGame()
      # Game logic.
      previous_health = self.health
      self.updateScorePerKillFromQueue()
      self.updateObjects(dt)
      if self.health > 0:
         self.spawnRock(dt)
         self.selectRock()
         self.destroyRock()
      # Game end.
      elif previous_health > 0:
         self.panel.showStartView()
      # Rendering.
      self.panel.deallocateAll()
      self.renderObjects()
      self.panel.noMoreAllocations()
      # Reset input variables.
      self.panel.mouse_pressed = False
      self.panel.start_pressed = False

   def updateScorePerKill(self, percent, message):
      self.score_per_kill = int(self.score_per_kill * (1.0 + percent / 100.0))
      self.panel.printMessage('[Score +{}%] {}'.format(percent, message))

   def updateScorePerKillFromQueue(self):
      if len(self.score_queue) > 0 and (self.current_time - self.score_queue[-1][0] >= 1.0):
         # Oldest element is at least one second old.
         hits = self.score_queue.pop()
         if len(hits) == 3:
            self.updateScorePerKill(50, '3 hits in one second')
         elif len(hits) == 4:
            self.updateScorePerKill(100, '4 hits in one second')
         elif len(hits) == 5:
            self.updateScorePerKill(200, '5 hits in one second')

   def updateHitQueue(self):
      # Hits within last second.
      hits = [self.current_time]
      for hit in self.hit_queue:
         if (hit is not None) and (self.current_time - hit <= 1.0):
            hits.append(hit)
      # If at least 3 hits within last second, then queue score update.
      # The reason for a queue is to avoid single hits contributing
      # to multiple score updates.
      if len(hits) >= 3:
         # If some of these hits are already queued, then replace the queued
         # score update, but only if the number of hits is greater.
         # If equal or less, then dismiss the current score update.
         if len(self.score_queue) > 0 and (self.score_queue[0][0] in hits):
            if len(hits) > len(self.score_queue[0]):
               self.score_queue[0] = hits
         # If none of these hits are queued, then add a new queued score update.
         else:
            self.score_queue.insert(0, hits)
      # Update hit queue.
      self.hit_queue[3] = self.hit_queue[2]
      self.hit_queue[2] = self.hit_queue[1]
      self.hit_queue[1] = self.hit_queue[0]
      self.hit_queue[0] = self.current_time

   def recordHit(self):
      self.updateHitQueue()
      # Update score per kill when a certain number of consecutive hits is made.
      self.consecutive_hits += 1
      if (self.consecutive_hits % 10) == 0:
         self.updateScorePerKill(self.consecutive_hits, '{} consecutive hits'.format(self.consecutive_hits))
      # Update score.
      self.score += self.score_per_kill
      self.panel.setScore(self.score)

   def recordMiss(self):
      self.consecutive_hits = 0

   def applyDamage(self, hit_position):
      x,y,z = hit_position
      dist_from_center = math.sqrt(x*x + y*y + z*z)
      damage_scale = max(1.0 - dist_from_center / (TARGET_AREA_SIZE/2), 0.0) # Range [0,1].
      damage = int(damage_scale * (DAMAGE_MAX - DAMAGE_MIN) + DAMAGE_MIN)
      self.health -= damage
      self.panel.setHealth(self.health)
      self.panel.shakeView(damage * (VIEW_SHAKE_SCALE / DAMAGE_MAX))

   def updateObjects(self, dt):
      # Update panel.
      self.panel.update(dt)
      # Update rocks.
      for rock in self.rocks:
         if rock.alive:
            rock.update(dt)
            if not rock.alive:
               # Crashed into the ground. Apply damage.
               self.applyDamage(rock.target)
         if not rock.alive:
            self.unused_rock = rock
      # Update booms.
      for boom in self.booms:
         if boom.alive:
            boom.update(dt)
         if not boom.alive:
            self.unused_boom = boom
      # Update beam.
      self.beam.update(dt)

   def spawnRock(self, dt):
      self.spawn_delay -= dt
      if self.spawn_delay <= 0:
         rate_ix = int((self.current_time - self.start_time) / SPAWN_RATE_DURATION)
         if rate_ix < len(SPAWN_RATES):
            self.spawn_delay = SPAWN_RATES[rate_ix]
         else:
            self.spawn_delay = SPAWN_RATES[-1]
         if self.unused_rock is None:
            self.unused_rock = Rock(self.sphere)
            self.rocks.append(self.unused_rock)
            self.rendered_objects.append(self.unused_rock)
         self.unused_rock.spawn()

   def selectRock(self):
      for rock in self.rocks:
         if rock.alive and ((self.selected_rock is None) or (self.selected_rock.position[2] > rock.position[2])):
            # Weak perspective projection.
            x = rock.position[0] / rock.position[2]
            y = rock.position[1] / rock.position[2]
            radius = ROCK_RADIUS / rock.position[2]
            # Calculate distance.
            dx = x - self.panel.mouse_position[0]
            dy = y - self.panel.mouse_position[1]
            dist = math.sqrt(dx*dx + dy*dy)
            if dist <= radius:
               self.selected_rock = rock

   def destroyRock(self):
      if self.panel.mouse_pressed:
         if self.selected_rock is not None:
            if self.unused_boom is None:
               self.unused_boom = Boom()
               self.booms.append(self.unused_boom)
               self.rendered_objects.append(self.unused_boom)
            self.unused_boom.spawn(self.selected_rock.position, self.selected_rock.velocity)
            self.beam.spawn(self.selected_rock.position)
            self.selected_rock.destroy()
            self.selected_rock = None
            self.recordHit()
         else:
            self.recordMiss()

   def renderObjects(self):
      if self.selected_rock is not None:
         self.selected_rock.select()
      self.rendered_objects.sort(key = lambda obj: obj.position[2])
      for obj in self.rendered_objects:
         if obj.alive:
            obj.renderInPanel(self.panel)
      if self.selected_rock is not None:
         self.selected_rock.deselect()

#===============================================================================
class Application:

   def __init__(self):
      self.wnd = tk.Tk()
      self.wnd.title('Orbital Target Practice')
      self.wnd.columnconfigure(0, weight = 1)
      self.wnd.rowconfigure(0, weight = 1)
      self.panel = ViewPanel(self.wnd, 0, 0, 800, 800)
      self.game = GameLogic(self.panel)
      self.wnd.after(FRAME_DELAY_MS, self.runFrame)

   def runFrame(self):
      self.game.runFrame()
      self.wnd.after(FRAME_DELAY_MS, self.runFrame)

   def run(self):
      self.wnd.mainloop()

if __name__ == '__main__':
   app = Application()
   app.run()
