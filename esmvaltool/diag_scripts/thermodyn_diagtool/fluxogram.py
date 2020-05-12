"""FLUX DIAGRAM PRODUCTION.

Created on Tue Jun 19 16:41:47 2018.

@author: Valerio2

Copyright 2018 Florian Ulrich Jehn

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from matplotlib import pyplot as plt


class Fluxogram():
    """The diagram flux module.

    A class to draw and maintain all fluxes and storages from a model or
    some similiar kind of thing to be drawn as a sequence of storages
    and fluxes.
    """

    def __init__(self, max_flux, max_storage, grid_size=20):
        """Initialize a fluxogram. must be called with.

        The arguments are:
            - max_flux: aximum flux of all fluxes; needed for scaling
            - max_storage: maximum storages of all storages; needed for scaling
            - grid_size:grid_size for drawing the fluxogram, determines how big
            everything is. Fluxes and storages  scaled accordingly
            - storages: all the storages the fluxogram has (usually empy to
            begin with)
            - fluxes: all the fluxes the fluxogram has (usually empty to begin
            with).
        """
        self.storages = []
        self.fluxes = []
        self.max_flux = max_flux
        self.max_storage = max_storage
        self.grid_size = grid_size

    def add_storage(self, name, amount, order, offset):
        """Add a storage to the storages of the fluxogram."""
        self.storages.append(
            Storage(name, self.grid_size, len(self.storages), amount, order,
                    offset))

    def add_flux(self, name, from_storage, to_storage, amount):
        """Add a flux to the fluxes of the fluxogram."""
        self.fluxes.append(
            Flux(name, self.grid_size, from_storage, to_storage, amount))

    def update_all_storages(self, amounts):
        """Update the amount of all storages."""
        for storage, amount in zip(self.storages, amounts):
            storage.update_storage(amount)

    def update_all_fluxes(self, amounts):
        """Update the amount of all fluxes."""
        for flux, amount in zip(self.fluxes, amounts):
            flux.update_flux(amount)

    def update_everything(self, amounts_storages, amounts_fluxes):
        """Update all fluxes and storages."""
        self.update_all_fluxes(amounts_fluxes)
        self.update_all_storages(amounts_storages)

    def draw(self, filen, listv):
        """Draw all fluxes and storages."""
        fig = plt.figure()
        frame1 = plt.axes()
        fig.set_size_inches(18.5, 10.5)
        # find the smallest/largest offset_ so the fluxogram can be drawn big
        # enough
        largest_offset = 0
        smallest_offset = 0
        largest_order = 0
        for storage in self.storages:
            if storage.offset > largest_offset:
                largest_offset = storage.offset
            if storage.offset < smallest_offset:
                smallest_offset = storage.offset
            if storage.order > largest_order:
                largest_order = storage.order
        # set y and x limits
        y_max = 0
        y_min = (largest_order + 1) * 2 * self.grid_size * -1
        x_max = (largest_offset + 2) * 2 * self.grid_size
        x_min = (smallest_offset - 1) * 2 * self.grid_size
        plt.axis([x_min, x_max, y_min, y_max])
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        # draw all fluxes
        dict_r = {
            'AZ+': listv[0],
            'ASE+': listv[2],
            'ATE+': listv[4],
            'A2KS': listv[6],
            'A2KT': listv[7],
            'KTE-': listv[8],
            'KSE-': listv[10],
            'KZ-': listv[12]
        }
        dict_oth = {
            'l': listv[14],
            'dn': listv[15],
            'rdn': listv[16],
            'ldn': listv[17],
            'up': listv[18],
            'lup': listv[19],
            'rup': listv[20]
        }
        switcher = {
            'l': self.leftarr_txt,
            'dn': self.dnarr_txt,
            'rdn': self.rdnarr_txt,
            'ldn': self.ldnarr_txt,
            'up': self.uparr_txt,
            'lup': self.luparr_txt,
            'rup': self.ruparr_txt
        }
        for flux in self.fluxes:
            idb = flux.name
            # scale the amount
            scaled_amount_flux = self.scaler(flux.amount, self.max_flux)
            # width multiplied  because if not, the arrows are so tiny
            arrow = plt.Arrow(
                flux.x_start,
                flux.y_start,
                flux.d_x,
                flux.d_y,
                width=scaled_amount_flux * 1.7,
                alpha=0.8)
            if flux.dire == 'r':
                for key in dict_r:
                    value = dict_r[key]
                    if idb == key:
                        plt.text(
                            flux.x_start + 0.25 * self.grid_size,
                            flux.y_start + 0.05 * self.grid_size,
                            value,
                            size=self.grid_size * 0.7)
            else:
                for key in dict_oth:
                    value = dict_oth[key]
                    if flux.dire == key:
                        switcher[flux.dire](value, flux, plt)
            plt.gca().add_patch(arrow)
        # draw all storages
        for storage in self.storages:
            # scale the amount
            scaled_amount_stor = self.scaler(storage.amount, self.max_storage)
            if scaled_amount_stor == 0:
                scaled_amount_stor = 0.0001
            # change_x and y, so the storages are centered to the middle
            # of their position and not to upper left
            x_p = (
                storage.x_p +
                (1 - storage.amount / self.max_storage) * 1.3 * self.grid_size)
            y_p = (
                storage.y_p -
                (1 - storage.amount / self.max_storage) * 1.3 * self.grid_size)
            rectangle = plt.Rectangle((x_p, y_p),
                                      scaled_amount_stor,
                                      -scaled_amount_stor,
                                      alpha=0.4)
            # label all storages
            plt.text(
                storage.x_p + 0.6 * self.grid_size,
                storage.y_p - 0.65 * self.grid_size,
                storage.name,
                fontsize=0.7 * self.grid_size)
            dict_s = {
                'AZ': listv[1],
                'ASE': listv[3],
                'ATE': listv[5],
                'KTE': listv[9],
                'KSE': listv[11],
                'KZ': listv[13]
            }
            for key in dict_s:
                value = dict_s[key]
                if storage.name == key:
                    plt.text(
                        storage.x_p + 0.6 * self.grid_size,
                        storage.y_p - 0.85 * self.grid_size,
                        value,
                        fontsize=0.7 * self.grid_size)
            # draw a date
            plt.gca().add_patch(rectangle)
        plt.savefig(filen)
        plt.close(fig)

    def dnarr_txt(self, value, flux, pltt):
        """Write text on arrow pointing down."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start - 0.2 * self.grid_size,
            y_start - 0.45 * self.grid_size,
            value,
            size=self.grid_size * 0.7,
            rotation=-90)

    def leftarr_txt(self, value, flux, pltt):
        """Write text on arrow pointing left."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start - 1.35 * self.grid_size,
            y_start + 0.05 * self.grid_size,
            value,
            size=self.grid_size * 0.7)

    def ldnarr_txt(self, value, flux, pltt):
        """Write text on arrow pointing down-left."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start - 0.35 * self.grid_size,
            y_start - 0.25 * self.grid_size,
            value,
            size=self.grid_size * 0.5,
            rotation=-110)

    def luparr_txt(self, value, flux, pltt):
        """Write text on arrow pointing up-left."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start - 0.35 * self.grid_size,
            y_start + 0.45 * self.grid_size,
            value,
            size=self.grid_size * 0.5,
            rotation=110)

    def rdnarr_txt(self, value, flux, pltt):
        """Write text on arrow pointing down-right."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start + 0.05 * self.grid_size,
            y_start - 0.25 * self.grid_size,
            value,
            size=self.grid_size * 0.5,
            rotation=-75)

    def ruparr_txt(self, value, flux, pltt):
        """Write text on arrow pointing up-right."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start - 0.1 * self.grid_size,
            y_start + 0.45 * self.grid_size,
            value,
            size=self.grid_size * 0.5,
            rotation=75)

    def uparr_txt(self, value, flux, pltt):
        """Write text on arrow pointing up."""
        x_start = flux.x_start
        y_start = flux.y_start
        pltt.text(
            x_start + 0.05 * self.grid_size,
            y_start + 0.75 * self.grid_size,
            value,
            size=self.grid_size * 0.7,
            rotation=90)

    def scaler(self, value_in, base_max):
        """Scale the values in the blocks of the diagram.

        Scale the fluxes and storages, so they don't overstep their
        grafical bounds must be called with:
            - valueIn: the value that needs rescaling
            - baseMax: the upper limit of the original dataset
                    ~ 100 for fluxes, ~250 for stores (in my model).
        """
        # baseMin: the lower limit of the original dataset (usually zero)
        base_min = 0
        # limitMin: the lower limit of the rescaled dataset (usually zero)
        limit_min = 0
        # limitMax: the upper limit of the rescaled dataset (in our case  grid)
        limit_max = self.grid_size
        # prevents wrong use of scaler
        if value_in > base_max:
            raise ValueError("Input value larger than base max")
        return (((limit_max - limit_min) * (value_in - base_min) /
                 (base_max - base_min)) + limit_min)


class Flux:
    """Contain a flux of a fluxogram."""

    def __init__(self, name, grid_size, from_storage, to_storage, amount=0):
        """Initialize a flux.

        Arguments are:
            - name: name of the flux
            - grid_size: grid size of the diagram
            - from_storage: storage the flux is originating from
            - to_storage: storage the flux is going into
            - amount: how much stuff fluxes.
        """
        self.name = name
        self.from_storage = from_storage
        self.to_storage = to_storage
        self.amount = amount
        self.grid_size = grid_size
        (self.x_start, self.y_start, self.x_end, self.y_end, self.d_x,
         self.d_y, self.dire) = (self.calc_start_end_dx_dy())

    def update_flux(self, amount):
        """Update the amount of the flux."""
        self.amount = amount

    def calc_start_end_dx_dy(self):
        """Scale the arrows.

        Calculate the starting and ending point of an arrow depending on the
        order and offset of the starting and ending storages. This helps
        determine the direction of the arrow
        returns the start and end xy coordinates of the arrow as tuples.
        """
        # arrow pointing to left up
        if (self.from_storage.offset > self.to_storage.offset
                and self.from_storage.order > self.to_storage.order):
            x_start = self.from_storage.x_p + 0.85 * self.grid_size
            y_start = self.from_storage.y_p - self.grid_size * 0.5
            x_end = self.to_storage.x_p + self.grid_size * 0.65
            y_end = self.to_storage.y_p - 0.7 * self.grid_size
            d_x = abs(x_start - x_end) * (-1)
            d_y = abs(y_start - y_end)
            dire = 'lup'
        # arrow pointing up
        elif (self.from_storage.offset == self.to_storage.offset
              and self.from_storage.order > self.to_storage.order):
            x_start = self.from_storage.x_p + 0.85 * self.grid_size
            y_start = self.from_storage.y_p - 0.5 * self.grid_size
            x_end = self.to_storage.x_p + 0.85 * self.grid_size
            y_end = self.to_storage.y_p - 0.25 * self.grid_size
            d_x = abs(x_start - x_end)
            d_y = abs(y_start - y_end)
            dire = 'up'
        # arrow pointing right up
        elif (self.from_storage.offset < self.to_storage.offset
              and self.from_storage.order > self.to_storage.order):
            x_start = (self.from_storage.x_p + self.grid_size)
            y_start = self.from_storage.y_p - 0.5 * self.grid_size
            x_end = self.to_storage.x_p + 0.05 * self.grid_size
            y_end = self.to_storage.y_p - 0.75 * self.grid_size
            d_x = abs(x_start - x_end)
            d_y = abs(y_start - y_end)
            dire = 'rup'
        # arrow pointing right
        elif (self.from_storage.offset < self.to_storage.offset
              and self.from_storage.order == self.to_storage.order):
            x_start = (self.from_storage.x_p + self.grid_size)
            y_start = self.from_storage.y_p - 0.8 * self.grid_size
            x_end = self.to_storage.x_p + 1.25 * self.grid_size
            y_end = self.to_storage.y_p - 0.8 * self.grid_size
            d_x = abs(x_start - x_end)
            d_y = abs(y_start - y_end)
            dire = 'r'
        # arrow pointing right down
        elif (self.from_storage.offset < self.to_storage.offset
              and self.from_storage.order < self.to_storage.order):
            x_start = (self.from_storage.x_p + 0.85 * self.grid_size)
            y_start = self.from_storage.y_p - 1.12 * self.grid_size
            x_end = self.to_storage.x_p + 0.85 * self.grid_size
            y_end = self.to_storage.y_p - 0.9 * self.grid_size
            d_x = abs(x_start - x_end)
            d_y = abs(y_start - y_end) * (-1)
            dire = 'rdn'
        # arrow pointing down
        elif (self.from_storage.offset == self.to_storage.offset
              and self.from_storage.order < self.to_storage.order):
            x_start = self.from_storage.x_p + 0.8 * self.grid_size
            y_start = (self.from_storage.y_p - 1.12 * self.grid_size)
            x_end = self.to_storage.x_p + 0.8 * self.grid_size
            y_end = self.to_storage.y_p - 1.4 * self.grid_size
            d_x = abs(x_start - x_end)
            d_y = abs(y_start - y_end) * (-1)
            dire = 'dn'
        # arrow pointing left down
        elif (self.from_storage.offset > self.to_storage.offset
              and self.from_storage.order < self.to_storage.order):
            x_start = self.from_storage.x_p + 0.75 * self.grid_size
            y_start = (self.from_storage.y_p - 1.1 * self.grid_size)
            x_end = self.to_storage.x_p + 0.6 * self.grid_size
            y_end = self.to_storage.y_p - 0.9 * self.grid_size
            d_x = abs(x_start - x_end) * (-1)
            d_y = abs(y_start - y_end) * (-1)
            dire = 'ldn'
        # arrow pointing left
        elif (self.from_storage.offset > self.to_storage.offset
              and self.from_storage.order == self.to_storage.order):
            x_start = self.from_storage.x_p + 0.5 * self.grid_size
            y_start = self.from_storage.y_p - 0.75 * self.grid_size
            x_end = self.to_storage.x_p + 0.25 * self.grid_size
            y_end = self.to_storage.y_p - 0.75 * self.grid_size
            d_x = abs(x_start - x_end) * (-1)
            d_y = abs(y_start - y_end)
            dire = 'l'
        # multiply by 0.9 so there is a gap between storages and arrows
        d_x = d_x * 0.75
        d_y = d_y * 0.75
        return x_start, y_start, x_end, y_end, d_x, d_y, dire


class Storage:
    """Contain a storage of a fluxogram."""

    def __init__(self, name, grid_size, number, amount=0, order=0, offset=0):
        """Initialize a storage.

        Arguments are:
                - name: name of the storage
                - number: consecutive number
                - grid_size of the diagram
                - amount: how much stuff is in it
                - order: how much down it is in the hierachie (starts with 0)
                - offset = how much the storage is offset to the left/right
                    in relationship to the center.
        """
        self.name = name
        self.amount = amount
        self.number = number
        self.order = order
        self.offset = offset
        self.grid_size = grid_size
        self.x_p, self.y_p = self.calculate_xy()

    def update_storage(self, amount):
        """Update the amount of the storage."""
        self.amount = amount

    def calculate_xy(self):
        """Provide coordinates of the blocks in the diagram.

        Calculate the xy coordinates of the starting point from where
        the rectangle is drawn. The additional multiplication by two is
        to produce the gaps in the diagram.
        """
        x_p = self.offset * self.grid_size * 2
        # multiply by -1 to draw the diagram from top to bottom
        y_p = self.order * self.grid_size * 2 * -1
        return x_p, y_p
