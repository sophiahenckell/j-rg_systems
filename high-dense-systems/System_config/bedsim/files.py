# Checked - A
from itertools import accumulate, chain, tee, product
from math import acos
from cmath import exp
from collections import deque

import numpy as np
import h5py

from bedsim.cell import DelaunayCellspace
import bedsim.boundary


class arb(object):
    pass


class ParticleFormat(object):
    def __init__(self, handle, num_timesteps):
        """
        @param handle Handle to the particle folder.
        @param tabsize Size of the table, usually #timesteps * #part. * #part properties
        """
        self.handle = handle
        self.num_timesteps = num_timesteps


class ShortFormat(ParticleFormat):
    def __flatten(self, list_of_lists):
        res = []
        for sublist in list_of_lists:
            if isinstance(sublist, list) or isinstance(sublist, np.ndarray):
                for val in sublist:
                    res.append(val)
            else:
                res.append(sublist)
        return res

    def get_number_of_timesteps(self):
        return self.get_number_of_frames()

    def get_number_of_frames(self):
        dynamics_handle = self.handle['/particles/dynamics']
        return dynamics_handle.shape[0]

    def load_statics(self):  # FIXME: load 'type'
        """ Set the data handle """
        statics_handle = self.handle['/particles/statics']

        """ Reconstruct column names and dimensions """
        column_names = statics_handle.attrs[
            'column_names']  # FIXME: these are in binary format, not uft-8! e.g. b'major' (could cause problems later)
        column_dimensions = statics_handle.attrs['column_dimensions']
        dim_indices = np.cumsum(column_dimensions)

        ids = statics_handle.attrs['particle_ids']
        types = statics_handle.attrs['particle_types']
        i2t = {}
        for i, t in zip(ids, types):
            # i2t[int(i)] = t.decode('utf8')
            i2t[int(i)] = t

        statics_bmf = []
        for static_dataset in statics_handle:  # each static_dataset represents the static values for one particle
            particles_dim_extract = np.split(static_dataset, dim_indices)[
                                    0:-1]  # omit last element which is for some reason an empty array
            particles_dim_extract_flattened = [x[0] if len(x) == 1 else x for x in
                                               particles_dim_extract]  # convert one dimensional arrays to numbers

            particle_dict = {}
            for k, v in zip(column_names, particles_dim_extract_flattened):
                # particle_dict[k.decode('utf8')] = v # construct the BMF dict for timestep 'step' for one particle
                particle_dict[k] = v  # construct the BMF dict for timestep 'step' for one particle

            # print("pdict id=",particle_dict['id'])
            ptype = i2t[particle_dict['id']]
            particle_dict['type'] = ptype

            statics_bmf.append(particle_dict)
            # yield particle_dict
        # print("sbmf=", statics_bmf)
        # return statics_bmf
        yield statics_bmf

    def load_dynamics(self):
        """
        Reads short format particle data and yields per time step in bedsim meta format (BMF).
        """

        """ Set the data handle """
        dynamics_handle = self.handle['/particles/dynamics']

        """ Reconstruct column names and dimensions """
        column_names = dynamics_handle.attrs[
            'column_names']  # FIXME: these are in binary format, not uft-8! e.g. b'major' (could cause problems later)
        column_dimensions = dynamics_handle.attrs['column_dimensions']
        dim_indices = np.cumsum(column_dimensions)  # later used for splitting

        """ get maximum time step and iterate over all times """
        times = dynamics_handle.shape[0]
        for step in range(times):  # FIXME: check if 'for data in dynamics_handle:' would also be possible
            data = dynamics_handle[step]  # get the raw data of timestep 'step'
            particles_bmf = []  # list of all particles' data for the timestep 'step'
            for particle_data in data:  # 'particle_data' is a list with one particle's properties
                particles_dim_extract = np.split(particle_data, dim_indices)[
                                        0:-1]  # omit last element which is for some reason an empty array
                particles_dim_extract_flattened = [x[0] if len(x) == 1 else x for x in
                                                   particles_dim_extract]  # convert one dimensional arrays to numbers
                particle_dict = {}
                for k, v in zip(column_names, particles_dim_extract_flattened):
                    particle_dict[k.decode('utf8')] = v  # construct the BMF dict for timestep 'step' for one particle
                particles_bmf.append(
                    particle_dict)  # append properties for this particle to the list of all particles' properties for that step
            yield particles_bmf
            # print("\n\n>>> BMF=", particles_bmf)

    def load_last_dynamics(self):
        """
        Reads short format particle data and returns last dynamics step in bedsim meta format (BMF).
        """
        if '/particles/dynamics' in self.handle:
            """ Set the data handle """
            dynamics_handle = self.handle['/particles/dynamics']
            if '/system/properties/summary_timestep_number' in self.handle:
                last_step = self.handle['/system/properties/summary_timestep_number'].value

                """ Reconstruct column names and dimensions """
                column_names = dynamics_handle.attrs[
                    'column_names']  # FIXME: these are in binary format, not uft-8! e.g. b'major' (could cause problems later)
                column_dimensions = dynamics_handle.attrs['column_dimensions']
                dim_indices = np.cumsum(column_dimensions)  # later used for splitting

                data = dynamics_handle[last_step]  # get the raw data of timestep 'step'
                particles_bmf = []  # list of all particles' data for the timestep 'step'
                for particle_data in data:  # 'particle_data' is a list with one particle's properties
                    particles_dim_extract = np.split(particle_data, dim_indices)[
                                            0:-1]  # omit last element which is for some reason an empty array
                    particles_dim_extract_flattened = [x[0] if len(x) == 1 else x for x in
                                                       particles_dim_extract]  # convert one dimensional arrays to numbers
                    particle_dict = {}
                    for k, v in zip(column_names, particles_dim_extract_flattened):
                        particle_dict[
                            k.decode('utf8')] = v  # construct the BMF dict for timestep 'step' for one particle
                    particles_bmf.append(
                        particle_dict)  # append properties for this particle to the list of all particles' properties for that step
                return particles_bmf
            else:
                return None

    def load_dynamics_step(self, i):
        """
        Reads short format particle data and returns i-th dynamics step in bedsim meta format (BMF).
        """
        if '/particles/dynamics' in self.handle:
            """ Set the data handle """
            dynamics_handle = self.handle['/particles/dynamics']
            if dynamics_handle[i]:
                """ Reconstruct column names and dimensions """
                column_names = dynamics_handle.attrs[
                    'column_names']  # FIXME: these are in binary format, not uft-8! e.g. b'major' (could cause problems later)
                column_dimensions = dynamics_handle.attrs['column_dimensions']
                dim_indices = np.cumsum(column_dimensions)  # later used for splitting

                data = dynamics_handle[i]  # get the raw data of timestep 'step'
                particles_bmf = []  # list of all particles' data for the timestep 'step'
                for particle_data in data:  # 'particle_data' is a list with one particle's properties
                    particles_dim_extract = np.split(particle_data, dim_indices)[
                                            0:-1]  # omit last element which is for some reason an empty array
                    particles_dim_extract_flattened = [x[0] if len(x) == 1 else x for x in
                                                       particles_dim_extract]  # convert one dimensional arrays to numbers
                    particle_dict = {}
                    for k, v in zip(column_names, particles_dim_extract_flattened):
                        particle_dict[
                            k.decode('utf8')] = v  # construct the BMF dict for timestep 'step' for one particle
                    particles_bmf.append(
                        particle_dict)  # append properties for this particle to the list of all particles' properties for that step
                return particles_bmf
            else:
                return None

    ### FIXME: integrate into save_to_file
    def save_statics(self, data):
        num_particles = len(data)
        num_properties = len(self.__flatten(data[0].values())) - 1  # subtract name

        attrs = {}
        arr = []
        for pdata in data:
            attrs[pdata['id']] = pdata.pop('type')
            arr.append(self.__flatten(pdata.values()))
            column_names = pdata.keys()  # ... overhead, I know, but only words that way, because python can't copy dicts without messing up their order...
            column_dimensions = [np.array(v).size for v in pdata.values()]
        static_data = np.array(arr)

        self.handle.require_group('/particles')  # ensure that the '/particles' group exists
        self.__dset_statics = self.handle['/particles'].require_dataset("statics", data=static_data,
                                                                        shape=(num_particles, num_properties),
                                                                        dtype='float64', compression="gzip",
                                                                        compression_opts=9)  # create dataset for static data

        self.__dset_statics.dims[0].label = 'particle'
        self.__dset_statics.dims[1].label = 'property'

        # NOTE: dict().copy() does not preserver dict order!!!
        # example_dict = data[0].copy() # get a copy of one entry to reconstruct the data structure properties
        # column_names = example_dict.keys()
        # column_dimensions = [np.array(v).size for v in example_dict.values()]

        self.__dset_statics.attrs['column_names'] = np.array([a.encode('utf8') for a in column_names],
                                                             dtype=h5py.special_dtype(
                                                                 vlen=str))  # save column names as attribute # FIXME: check if sorting is correct
        self.__dset_statics.attrs['column_dimensions'] = column_dimensions
        # print("cn=",column_names)
        # print("cd=",column_dimensions)

        # own attribute for each particle.... messy...
        # for k,v in attrs.items():
        #    self.__dset_statics.attrs[str(k)] = v

        # save attrs as string list dataset... not working atm... this is messy in h5   
        # dt = h5py.special_dtype(vlen=str)
        # type_handle = self.handle['/particles'].require_dataset('types', data=np.asarray( [x.encode('utf8') for x in list(attrs.values())] ), shape=(1,len(attrs)), dtype=dt)

        # save as attribute... might be bigger than 64k but should work (c.f. https://www.hdfgroup.org/HDF5/doc/UG/13_Attributes.html)
        ids = attrs.keys()
        types = attrs.values()
        self.__dset_statics.attrs['particle_ids'] = np.array([str(id).encode('utf8') for id in ids],
                                                             dtype=h5py.special_dtype(vlen=str))
        self.__dset_statics.attrs['particle_types'] = np.array([a.encode('utf8') for a in types],
                                                               dtype=h5py.special_dtype(vlen=str))

    def save_dynamics(self, data, time):
        """
        NOTE: if 'time' if not an integer this will not work!!
        """
        if self.__dset is None:
            num_particles = len(data)
            # print("\n\n<<< PROPS: ", self.__flatten(data[0].values()))
            num_properties = len(self.__flatten(data[0].values())) - 1  # subtract 'type'
            self.handle.require_group('/particles')
            self.__dset = self.handle['/particles'].require_dataset("dynamics", shape=(
            int(self.num_timesteps), num_particles, num_properties), dtype='float64', compression="gzip",
                                                                    compression_opts=9)

            self.__dset.dims[0].label = 'time'
            self.__dset.dims[1].label = 'particle'
            self.__dset.dims[2].label = 'property'

            """ Save column names as attribute """
            # FIXME: ensure that dict sorting is stable!!! => seems like it...
            """
            example_dict = data[0].copy() # get a copy of one entry to reconstruct the data structure properties
            example_dict.pop('type')
            column_names = example_dict.keys()
            column_dimensions = [np.array(v).size for v in example_dict.values()]
            """
            column_names = list(data[0].keys())
            typepos = column_names.index(
                'type')  # get the position of the type property to remove it from the dimensions list
            column_names.remove('type')
            column_dimensions = [np.array(v).size for v in data[0].values()]
            del column_dimensions[typepos]

            self.__dset.attrs['column_names'] = [a.encode('utf8') for a in column_names]
            self.__dset.attrs['column_dimensions'] = column_dimensions

        arr = []
        for pdata in data:  ### FIXME: check what is more efficient: using id as dimenaion or table index
            pdata.pop('type')  # particle type is not saved each time
            """ NOTE: dict.values() returns a dict_keys([..]) object => iterable like a list but not accessible like a list (like e.g. l[1]) """
            arr.append(self.__flatten(pdata.values()))

            # self.__dset[time,pdata['id']] = self.__flatten(pdata.values())
            # self.__dset[pdata['id'], time] = self.__flatten(pdata.values())

        # print("arr=", np.array(arr).shape)
        self.__dset[time] = np.array(arr)

    def __init__(self, handle, num_timesteps):
        ParticleFormat.__init__(self, handle, num_timesteps)
        self.__dset = None
        self.__dset_statics = None


class LongFormat(ParticleFormat):  # aka verbose format
    """
    The long format class defines a very verbose saving schema.
    """

    def get_number_of_timesteps(self):
        return self.get_number_of_frames()

    def get_number_of_frames(self):
        handle_time_folder = self.handle['/particles/dynamics']
        time_list = list(handle_time_folder.keys())
        return len(time_list)

    #### FIXME: integrate into the standard loading method
    #### THIS IS MORE OR LESS A COPY OF THE FUNCTION BELOW
    def load_statics(self):
        handle_folder = self.handle['/particles/statics']
        particle_list = list(handle_folder.keys())
        particles = []
        for gid in particle_list:
            (type, id) = gid.split('_')
            phandle = handle_folder[gid]
            """ Convert the HDF dataset to python/numpy types """
            # TODO: Questions:
            #  (1) is this really necessary and
            #  (2) can this be done any faster?
            #  (3) The problem is that phandle.items() is in a special h5py format... most likely not directly compatible to the types expected by the rest of the app
            pppdata = {'id': id, 'type': type}
            for (k, v) in phandle.items():
                pppdata[k] = v.value
            particles.append(pppdata)
        yield particles

    def load_dynamics(self):  # FIXME!
        """
        Returns an iterator over particle times providing particle data.
        """
        # FIXME: extract type from gid
        handle_time_folder = self.handle['/particles/dynamics']
        time_list = list(handle_time_folder.keys())
        time_list = sorted([float(time) for time in
                            time_list])  # assure that time is ascendant and not influenced by string ordering, e.g. '10.0' < '2.0'

        for time in time_list:
            if time == 0.0:  # we used an int to save 0 state and converted to float above
                time = '0'  # FIXME: this can be done in a nicer way i guess :P
            handle_particle_folder = self.handle['/particles/dynamics/%s' % time]
            particle_list = list(handle_particle_folder.keys())
            particles = []
            for gid in particle_list:
                # CHECK FOR PARTICLE TYPE!
                (type, id) = gid.split('_')
                phandle = handle_particle_folder[gid]
                """ Convert the HDF dataset to python/numpy types """
                # TODO: Questions:
                #  (1) is this really necessary and
                #  (2) can this be done any faster?
                #  (3) The problem is that phandle.items() is in a special h5py format... most likely not directly compatible to the types expected by the rest of the app
                pppdata = {'id': id, 'type': type}
                for (k, v) in phandle.items():
                    pppdata[k] = v.value
                particles.append(pppdata)
            yield particles

    def load_dynamics_step(self, i):
        """
        Returns last dynamics.
        """
        # FIXME: extract type from gid
        handle_time_folder = self.handle['/particles/dynamics']

        handle_particle_folder = self.handle['/particles/dynamics/%s' % i]
        particle_list = list(handle_particle_folder.keys())
        particles = []
        for gid in particle_list:
            (type, id) = gid.split('_')
            phandle = handle_particle_folder[gid]
            pppdata = {'id': id, 'type': type}
            for (k, v) in phandle.items():
                pppdata[k] = v.value
            particles.append(pppdata)
        return particles

    def load_last_dynamics(self):
        """
        Returns last dynamics.
        """
        # FIXME: extract type from gid
        handle_time_folder = self.handle['/particles/dynamics']
        time_list = list(handle_time_folder.keys())
        time_list = sorted([float(time) for time in
                            time_list])  # assure that time is ascendant and not influenced by string ordering, e.g. '10.0' < '2.0'
        last = max(time_list)

        handle_particle_folder = self.handle['/particles/dynamics/%s' % last]
        particle_list = list(handle_particle_folder.keys())
        particles = []
        for gid in particle_list:
            (type, id) = gid.split('_')
            phandle = handle_particle_folder[gid]
            pppdata = {'id': id, 'type': type}
            for (k, v) in phandle.items():
                pppdata[k] = v.value
            particles.append(pppdata)
        return particles

    def save_statics(self, data):
        self.save_dynamics(data, None)

    def save_dynamics(self, data, time=None):  ### FIXME:...
        """
        @param time Time step to save data for. If time is 'None' data is written as static data
        @param data Data should be of the form [{'name':'Particle_01', 'position':array([0,0]), ...}]
        """
        if time is not None:  ## FIXME: maybe there is a nicer way to do this... like kwargs or something
            self.handle.require_group('/particles/dynamics/%s' % time)
            thandle = self.handle['/particles/dynamics/%s' % time]
        else:
            self.handle.require_group('/particles/statics')
            thandle = self.handle['/particles/statics']
        for pdata in data:
            # group_name = pdata.pop('name')
            group_name = "%s_%s" % (pdata.pop('type'), pdata.pop('id'))
            thandle.require_group(group_name)
            phandle = thandle[group_name]

            for key, value in pdata.items():
                shape = () if np.isscalar(value) else np.array(value).shape
                try:
                    phandle.require_dataset('%s' % key, data=value, shape=shape, dtype='float64')
                except RuntimeError:
                    print("Saving runtime error: maybe duplicate key: ", key, " at time ", time, ".")


class Bedfile(object):
    """
    Brownian event driven simulations IO library.
    Manages saving and loading configuration and simulation files.
    
    A bedfile is a hdf5 file with the following structure:
    root
     |-- system
     |    |-- boundary
     |    +-- cellspace
     |    |    |-- grid
     |    +-- properties
     |         |-- localtime
     +-- particles
          |-- dynamics
          |-- statics
     
     
    Inter-Module particle data exchange is performed in a meta format (BMF):
     
    > For a certain t: [{'property0':value0, 'prop1':val1, ...}, ...]
    
    The outer list carries dictionaries for each particle. Each dictionary carries all
    the information about the particles (name, position, ...).
    """

    def __init__(self, filename):
        self.handle = h5py.File(filename, "a")
        self.system = BedfileSystem(self.handle)


class BedfileReader(Bedfile):
    def _get_format(self):
        print("inside BedfileReader, getting format")
        return self.handle['/particles'].attrs['format']  # .value.decode('utf-8')

    # FIXME: this function is just an ugly hack for the 'antrag' series... improve it in the future
    def calc_packing_fraction(
            self):  # FIXME: variable area # FIXME: THIS IS A QUICK AND DIRTY IMPLEMENTATION FOR THE ANTRAG SERIES...
        # FIXME: box should be better because P3 is partly radial symmetric
        epsilon = 0.5
        packing_fraction = []  # packing fraction for each time step (because we look at a certain interval only)
        (xmin, xmax, ymin, ymax, zmin, zmax) = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
        center = np.array([(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2])
        # area = (xmax - xmin) * (ymax - ymin) # FIXME! # box area
        volume = 4. * np.pi * epsilon ** 3 / 3.  # circular area
        data = self.particles.load_dynamics()  # FIXME: geometric properties will be stored in statics table! Adjust this!
        for t in data:
            particle_volume = 0
            for particle in t:
                if np.linalg.norm(particle['position'] - center) < epsilon:  # circle
                    particle_volume += 4. * np.pi * particle['major'] * particle['minor'] * particle['minor2'] / 3.
            packing_fraction.append(particle_volume / volume)
        return packing_fraction



    def __pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)

    def msd(self):
        """
        Calculate MSD and revert boundary conditions or angular mod 2pi
        """
        print("inside msd")
        arbsys = arb()  # set up a dummy system... needed to use boundary module... maybe fix this...
        arbsys.cellspace = DelaunayCellspace(system=arbsys,
                                             grid_points=self.system.grid)  ## this is needed to use boundary... not nice
        boundary = getattr(bedsim.boundary, self.system.boundary)(
            system=arbsys)  # now boundary.delta_dir(r1,r2) can be used

        data = self.particles.load_dynamics()  # list of times with list of particles with list of properties
        a = next(
            data)  # PROBLEM: at the moment we write pos=0 for t=0 in dynamics which is not good. # FIXME: once this is fixed we don't need to concat the statics data anymore!!!! #######

        d0 = next(self.particles.load_statics())
        data = [d0] + list(
            data)  # this is needed, because pos=0 for t=0 in dynamics data! ##### FIXME: REMOVE ONCE DYNAMICS SAVES POSITION PROPERLY FOR t=0!! #####

        displacement = []  # i.e. cumulative positions
        for t1, t2 in self.__pairwise(data):

            if not t1[0]['pinned'] or t1[0]['pinned'] == 0:
                disp_per_time = []

                for t1p, t2p in zip(t1, t2):  # iterate over each particle of the two timesteps combined
                    disp_per_time.append(
                        (t1p['cumulative_position'], t2p['cumulative_position']))  # calc shift during last timestep
                if displacement != []:
                    dp2 = []
                    # now calculate the cumulativep position by adding disp_per_time (last shift) to last cumulative position (displacement[-1])
                    for i, j in zip(displacement[-1], disp_per_time):
                        dp2.append(i + j)
                    displacement.append(dp2)  # finally append the cumulative value
                else:
                    displacement.append(disp_per_time)

        # calc msd
        msd = []
        for t in displacement:
            msd_t = 0
            for p_traj in t:
                msd_t += np.dot(p_traj, p_traj)
            msd.append(msd_t / len(t))
        return msd

    #     def msd_precalc(self):
    #         """
    #         MSD calcumation from pre-calculated cumulative positions
    #         """
    #         #static_data = self.particles.load_statics()
    #         data = self.particles.load_dynamics() # list of times with list of particles with list of properties
    #
    #         a = next(data)
    #
    #         msd = []
    #         for particles_at_t in data:
    #             csum2 = 0
    #             unpinned_particles = 0
    #             for particle_t in particles_at_t:
    #                 if not particle_t['pinned']:
    #                     unpinned_particles += 1
    #                     cpos = particle_t['cumulative_position']
    #                     #print("cpos=",np.dot(cpos, cpos))
    #                     csum2 += np.dot(cpos, cpos)
    #             msd_t = csum2/unpinned_particles
    #             msd.append(msd_t)
    #         return msd





    def msr(self):
        print("inside msr")
        data = self.particles.load_dynamics()
        a = next(data)

        time = 0
        displacement = []
        for t1, t2 in self.__pairwise(data):
            if not t1[0]['pinned'] or t1[0]['pinned'] == 0:
                disp_per_time = []
                for t1p, t2p in zip(t1, t2):
                    theta1 = t1p['angle']
                    theta2 = t2p['angle']

                    z = acos(np.dot(theta2[:3], theta1[:3]))  # get only the vector part of the quaternion

                    disp_per_time.append(z)

                if displacement != []:
                    dp2 = []
                    for i, j in zip(displacement[-1], disp_per_time):
                        dp2.append(i + j)
                    displacement.append(dp2)
                else:
                    displacement.append(disp_per_time)
                time += 1

        # calc msr
        msr = []
        for t in displacement:
            msr_t = 0
            for ang_traj in t:
                msr_t += np.dot(ang_traj, ang_traj)  # ang_traj * ang_traj
            msr.append(msr_t / len(t))
        return msr

    def sk_static(self):
        """
        Calculate static distribution function
        """
        arbsys = arb()  # set up a dummy system... needed to use boundary module... maybe fix this...
        arbsys.cellspace = DelaunayCellspace(system=arbsys,
                                             grid_points=self.system.grid)  ## this is needed to use boundary... not nice
        boundary = getattr(bedsim.boundary, self.system.boundary)(
            system=arbsys)  # now boundary.delta_dir(r1,r2) can be used

        data = self.particles.load_dynamics()  # list of times with list of particles with list of properties
        #         a = next(data) # PROBLEM: at the moment we write pos=0 for t=0 in dynamics which is not good. # FIXME: once this is fixed we don't need to concat the statics data anymore!!!! #######
        data = list(data)
        d0 = data[len(data) - 1]
        #         d0 = next(self.particles.load_statics())
        #         data = [d0]+list(data) # this is needed, because pos=0 for t=0 in dynamics data! ##### FIXME: REMOVE ONCE DYNAMICS SAVES POSITION PROPERLY FOR t=0!! #####
        #         nn = len(data)
        #         N = len(data[0]) # number of particles in the system
        L = 20.92271492899523  # 7.846018098373213
        #         particle_pairs = product(list(range(N),repeat=2)

        q_space = np.linspace(start=0, stop=40, num=41,
                              endpoint=True).tolist()  # num+1 because endpoint counts as well...
        q_data = [[float(qx), float(qy), float(qz)] for qx in q_space for qy in q_space for qz in q_space]
        q_data = np.delete(q_data, 0, 0)
        sk_static = []

        rmin = []
        for particle in d0:  # get the statistics of the first one
            if particle['pinned'] == 0:
                rmin.append(boundary.unwrap(particle['position']))
                # start here

            #         seen = set()
            #         seenlist = []
            #         qlist = []
            #
            #         for (qx,qy,qz) in q_data:
            #             normq = np.linalg.norm([qx,qy,qz])
            #             if normq not in seen:
            #                 seen.add(normq)
            #         for (qx,qy,qz) in q_data:
            #             normq = np.linalg.norm([qx,qy,qz])
            #             ind = seenlist.index(normq)
            #             qlist[ind].append([qx,qy,qz])
            #
            #         qlist = deque(qlist)
            #         sk_static = []
            #         sk = []
            #
            #         for sets in qlist:
            #             sets = deque(sets)
            #             mag_q = sets.popleft() # normalized q
            #             rhoq2_per_time = 0.
            #             qnum = len(sets)
            #             while sets: # the remainder is a list of wave vectors
            #                 q = sets.popleft()
            # #                 rhoq_per_time = (np.sum( [ exp(-1j * np.dot(q,r) * 2*np.pi/L) for r in rmin ] ) )
            #                 cosq = np.sum( np.cos(np.dot(q,r) * 2.*np.pi/L) for r in rmin )
            #                 sinq = np.sum( np.cos(np.dot(q,r) * 2.*np.pi/L) for r in rmin )
            #                 rhoq_per_time = cosq*cosq + sinq*sinq
            #                 rhoq2_per_time += rhoq_per_time
            # #                 rhoq2_per_time += np.real(rhoq_per_time*np.conj(rhoq_per_time))/len(d0)
            #             sk_static.append([mag_q,rhoq2_per_time/qnum])
            #
            #         return sk_static

        qlist = np.zeros((int(np.sqrt(3. * len(q_space) ** 2)), 2))
        for q in q_data:
            normq = np.linalg.norm(q)

            if normq != 0:

                if normq % 1 == 0:
                    bin = qlist[int(normq)][0]
                    rhoq2_per_time = qlist[int(normq)][1]

                    bin += 1

                    rhoq_per_time = (np.sum([exp(-1j * np.dot(q, r) * 2 * np.pi / L) for r in rmin]))
                    rhoq2_per_time += np.real(rhoq_per_time * np.conj(rhoq_per_time)) / len(d0)

                    qlist[int(normq)][0] = bin
                    qlist[int(normq)][1] = rhoq2_per_time

        sk_static = []

        for (bins, rhos) in qlist.tolist():
            if bins != 0:
                sk_static.append(rhos / bins)
            else:
                sk_static.append(rhos)

        return sk_static

    #     def sk(self):
    #         """
    #         Calculate MSD and revert boundary conditions or angular mod 2pi
    #         """
    #         print("inside msd")
    #         arbsys = arb() # set up a dummy system... needed to use boundary module... maybe fix this...
    #         arbsys.cellspace = DelaunayCellspace(system=arbsys, grid_points=self.system.grid) ## this is needed to use boundary... not nice
    #         boundary = getattr(bedsim.boundary, self.system.boundary)(system=arbsys) # now boundary.delta_dir(r1,r2) can be used
    #
    #         data = self.particles.load_dynamics() # list of times with list of particles with list of properties
    #         a = next(data) # PROBLEM: at the moment we write pos=0 for t=0 in dynamics which is not good. # FIXME: once this is fixed we don't need to concat the statics data anymore!!!! #######
    #
    #         d0 = next(self.particles.load_statics())
    #         data = [d0]+list(data) # this is needed, because pos=0 for t=0 in dynamics data! ##### FIXME: REMOVE ONCE DYNAMICS SAVES POSITION PROPERLY FOR t=0!! #####
    #
    #         N = len(data[0]) # number of particles in the system
    #         box_size = 7.846018098373213
    #         particle_pairs = product(list(range(N),repeat=2)
    #
    #
    #
    #         q_space = np.linspace(start=0, stop=10, num=100, endpoint=True).tolist() # num+1 because endpoint counts as well...
    #         q_data = [ [float(qx),float(qy),float(qz)] for qx in q_space for qy in q_space for qz in q_space]
    #
    #
    #         displacement = [] # i.e. cumulative positions
    #         coherent = [] # coherent part of the structure factor
    #         incoherent = [] # incoherent part of the structure factor
    #         sum_coherent = 0.0
    #         sum_incoherent = 0.0
    #
    #         for t1,t2 in self.__pairwise(data):
    #
    #             if not t1[0]['pinned'] or t1[0]['pinned'] == 0:
    #                 disp_per_time = []
    #                 incoherent_per_time = []
    #                 coherent_per_time = []
    #
    #                 for (parti,partj) in particle_pairs:
    #
    #                     if parti == partj:
    #                         rti = boundary.delta_dir(t1[parti]['position'] - t2[partj]['position']) # calc shift during last timestep
    #                         incoherent_per_time.append(rti)
    #                     else:
    #                         rtc = boundary.delta_dir(t1[parti]['position'] - t2[partj]['position']) # calc shift during last timestep
    #                         coherent_per_time.append(rtc)
    #
    #                 f_coherent =  1/(2 * np.pi) * ( np.sum( [exp(+1j * np.dot(q,r) ) for q in q_data for r in coherent_per_time] ) )
    #                 f_incoherent =  1/(2 * np.pi) * ( np.sum( [exp(-1j * np.dot(q,r) ) for q in q_data for r in incoherent_per_time] ) )
    #
    # #                 if coherent_per_time != []:
    # #                     coh2 = []
    # #                     # now calculate the cumulativep position by adding disp_per_time (last shift) to last cumulative position (displacement[-1])
    # #                     for i,j in zip(coherent[-1], f_coherent):
    # #                         coh2.append(i+j)
    # #                     coherent.append(coh2) # finally append the cumulative value
    # #                 else:
    # #                     coherent.append(f_coherent)
    #
    #                 if incoherent_per_time != []:
    #                     inc2 = []
    #                     # now calculate the cumulativep position by adding disp_per_time (last shift) to last cumulative position (displacement[-1])
    #                     for i,j in zip(incoherent[-1], f_incoherent):
    #                         inc2.append(i+j)
    #                     incoherent.append(inc2) # finally append the cumulative value
    #                 else:
    #                     incoherent.append(f_incoherent)
    #
    #         # calc
    #         msd = []
    #         for t in incoherent:
    #             msd_t = 0
    #             for p_traj in t:
    #                 msd_t += np.dot(p_traj, p_traj)
    #             msd.append(msd_t/len(t))
    #         return msd


    '''
    def msr_precalc(self):
        """
        MSD implementation which uses precalculated cumulative angles
        """
        data = self.particles.load_dynamics()
        a = next(data)
        
        time = 0
        displacement = []
        for t1, t2 in self.__pairwise(data):
            if not t1[0]['pinned'] or t1[0]['pinned'] == 0:
                disp_per_time = []
                for t1p, t2p in zip(t1, t2):
                    # like delta_dir but independent from boundary...
                    angular_difference = t2p['cumulative_angle'] - t1p['cumulative_angle'] # FIXME: check order!
                    disp_per_time.append( angular_difference )

                if displacement != []:
                    dp2 = []
                    for i,j in zip(displacement[-1], disp_per_time):
                        dp2.append(i+j)
                    displacement.append(dp2)
                else:
                    displacement.append(disp_per_time)
                time += 1

        # calc msr
        msr = []
        for t in displacement:
            msr_t = 0
            for ang_traj in t:
                msr_t += ang_traj * ang_traj
            msr.append(msr_t/len(t))
        return msr
    '''

    #     def msr_precalc(self):
    #         static_data = self.particles.load_statics()
    #         data = self.particles.load_dynamics() # list of times with list of particles with list of properties
    #
    #         a = next(data)
    #
    #         msr = []
    #         for particles_at_t in data:
    #             msr_t = 0
    #             unpinned_particles = 0
    #             for particle_t in particles_at_t:
    #                 if not particle_t['pinned']:
    #                     unpinned_particles += 1
    #                     cang = particle_t['cumulative_angle']
    #                     msr_t += np.dot(cang, cang)
    #             msr.append(msr_t/unpinned_particles)
    #         return msr # FIXME: CONSIDER MOMENT OF INERTIA




    def ft(self):
        """
        Returns fourier transform of all simulation steps.
        The fourier transform is performed over delta functions at the particles' current position.
        """
        tdata = self.particles.load_dynamics()  # list of times with list of particles with list of properties
        a = next(tdata)  # omit first step

        for data in tdata:
            z = []
            for p in data:
                # z.append(p['cumulative_position'])
                if not p['pinned']:
                    z.append(p['position'])

            from cmath import exp
            z = np.array(z)
            ft = lambda w1, w2, w3: 1 / (2 * np.pi) * (
            np.sum([exp(-1j * (xyz[0] * w1 + xyz[1] * w2 + xyz[2] * w3)) for xyz in z]))
            yield ft


        #     def ft_static(self):
        #         """
        #         Returns fourier transform of a certain simulation step.
        #         The fourier transform is performed over delta functions at the particles' current position.
        #         """
        #         static_data = self.particles.load_statics()
        #         static_data = next(static_data) # FIXME: fix this...
        #         #data = self.particles.load_dynamics() # list of times with list of particles with list of properties
        #         #data = self.particles.load_last_dynamics()
        #
        #         #a = next(data)
        #
        #         z = []
        #         for p in static_data:
        #             #z.append(p['cumulative_position'])
        #             if not p['pinned']:
        #                 z.append(p['position'])
        #
        #         from cmath import exp
        #         z = np.array(z)
        #         ft = lambda w1,w2: 1/(2 * np.pi) * ( np.sum( [exp(-1j * (xy[0] * w1 + xy[1] * w2) ) for xy in z] ) )
        #
        #         grid_data = self.system.grid
        #         [x,y] = grid_data.transpose()
        #         [xmin, xmax, ymin, ymax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y)]
        #         [width, height] = [xmax-xmin, ymax-ymin]
        #         heaviside = lambda x: 1 if x>=0 else 0
        #         s1=width/2
        #         s2=height/2
        #         box_ft = lambda w1,w2: (2*heaviside(2*s1)*heaviside(2*s2)*np.sin(s1*w1)*np.sin(s2*w2))/(np.pi*w1*w2)
        #         diff_ft = lambda w1,w2: ft(w1,w2)/box_ft(w1,w2)
        #
        #
        #         return ft
        #         #return diff_ft


        #     ## FIXME: REMOVE THIS METHOD AND GENERALIZE THE METHOD ABOVE!!!
        #     def ft_last(self):
        #         data = self.particles.load_last_dynamics()
        #         z = []
        #         for p in data:
        #             if not p['pinned']:
        #                 z.append(p['position'])
        #         from cmath import exp
        #         z = np.array(z)
        #         ft = lambda w1,w2: 1/(2 * np.pi) * ( np.sum( [exp(-1j * (xy[0] * w1 + xy[1] * w2) ) for xy in z] ) )
        #
        #
        #         grid_data = self.system.grid
        #         [x,y] = grid_data.transpose()
        #         [xmin, xmax, ymin, ymax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y)]
        #         [width, height] = [xmax-xmin, ymax-ymin]
        #         heaviside = lambda x: 1 if x>=0 else 0
        #         s1=width/2
        #         s2=height/2
        #         box_ft = lambda w1,w2: (2*heaviside(2*s1)*heaviside(2*s2)*np.sin(s1*w1)*np.sin(s2*w2))/(np.pi*w1*w2)
        #         diff_ft = lambda w1,w2: ft(w1,w2)/box_ft(w1,w2)
        #
        #         return ft
        #         #return diff_ft


        #     def ft_step(self, i):
        #         # ft at step i
        #         data = self.particles.load_dynamics_step()
        #         z = []
        #         for p in data:
        #             if not p['pinned']:
        #                 z.append(p['position'])
        #         from cmath import exp
        #         z = np.array(z)
        #         ft = lambda w1,w2: 1/(2 * np.pi) * ( np.sum( [exp(-1j * (xy[0] * w1 + xy[1] * w2) ) for xy in z] ) )
        #         return ft

    def __init__(self, filename):
        Bedfile.__init__(self, filename)
        self.handle = h5py.File(filename, "r")
        fmt = self._get_format()
        self.particles = getattr(bedsim.files, fmt)(self.handle, None)


class BedfileWriter(Bedfile):
    def __init__(self, filename, particle_format="LongFormat", num_timesteps=None):
        """
        @param filename Path to the h5 file to write to.
        @param particle_format The format particles should be written in. 'LongFormat' is default for maximum compatibility. All children of the 'ParticleFormat' class are valid.
        @param num_timesteps Number of time steps of the whole simulation (used to reserve space). Omit if just saving static data. 
        """
        Bedfile.__init__(self, filename)
        self.handle = h5py.File(filename, "a")
        self.particles = getattr(bedsim.files, particle_format)(self.handle, num_timesteps)
        self.handle.require_group('/particles')
        self.handle['/particles'].attrs['format'] = particle_format


class BedfileSystem(object):
    @property
    def boundary(self):
        return self.handle['/system/boundary'].value  # .decode('utf-8')

    @boundary.setter
    def boundary(self, boundary_data):
        self.handle.require_group('/system')
        # self.handle.require_dataset('/system/boundary', data=np.string_(boundary_data), shape=(), dtype="S10")
        self.handle.require_dataset('/system/boundary', data=boundary_data, shape=(),
                                    dtype=h5py.special_dtype(vlen=str))

    @property
    def grid(self):
        return self.handle['/system/cellspace/grid'].value  # return grid points

    @grid.setter
    def grid(self, grid_data):
        self.handle.require_group('/system/cellspace')
        self.handle.require_dataset('/system/cellspace/grid', data=grid_data, shape=np.array(grid_data).shape,
                                    dtype="float64")

    @property
    def system_properties(self):
        properties_handle = self.handle['/system/properties']
        properties = properties_handle.items()
        # convert property hdf format to standard types # FIXME: is this necessary?
        pppdata = {}
        for (k, v) in properties:
            pppdata[k] = v.value
        return pppdata

    @system_properties.setter
    def system_properties(self, property_data):
        self.handle.require_group('/system/properties')
        """
        for k in property_data:
            #if property_data[k] is not None:
            value = property_data[k][0]
            shape = () if np.isscalar(value) or value is None else value.shape
            self.handle.require_dataset('/system/properties/%s' % k, data=value, shape=shape, dtype=property_data[k][1])
        """
        for k, v in property_data.items():
            if v is not None:
                shape = () if np.isscalar(v) or v is None else v.shape
                if '/system/properties/' + k in self.handle:
                    del self.handle['/system/properties/' + k]
                self.handle.require_dataset('/system/properties/%s' % k, data=v, shape=shape, dtype=type(v))

    @property
    def system(self):
        boundary = {'boundary': self.boundary()}
        grid = {'grid': self.grid_data()}
        properties = {'properties': self.system_properties()}
        return properties.update(boundary, grid)

    @system.setter
    def system(self, system_data):
        self.boundary(system_data['boundary'])
        self.grid(system_data['grid'])
        self.system_properties(system_data['properties'])

    def __init__(self, handle):
        self.handle = handle
