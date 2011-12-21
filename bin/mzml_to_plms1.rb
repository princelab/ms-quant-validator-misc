#!/usr/bin/env ruby

require 'trollop'
require 'mzml'
require 'ms/msrun/plms1'

parser = Trollop::Parser.new do
  banner = %Q{usage: #{File.basename(__FILE__)} <file>.mzML
output: <file>.plms1
}
  opt :filter, "accepts only non-overlapping rt filters of the form: <rt>:<rt>,<m/z>:<m/z>[,...]", :type => :string, :multi => true
end

opt = parser.parse(ARGV)

if ARGV.size == 0
  parser.educate && exit
end

class Filter
  attr_accessor :rt_range
  attr_accessor :mz_ranges
  def initialize(rt_range, *mz_ranges)
    @rt_range = rt_range
    @mz_ranges = mz_ranges
  end

  # <rt>:<rt>,<m/z>:<m/z>[,...]
  def self.from_string(string)
    ranges = string.split(',').map {|range| Range.new(*range.split(':').map(&:to_f)) }
    Filter.new ranges.shift, *ranges
  end

  # if it passes the rt range, will return an array of arrays (mz and
  # intensity) of any peaks passing the m/z filters.  returns nil if it
  # doesn't pass the rt range. ar_of_arrays is an array of two arrays, m/z's
  # and intensities.
  def filter(time, ar_of_arrays)
    if @rt_range === time
      #puts "ACCEPTING #{@rt_range} === #{time}"
      new_mzs = []
      new_int = []
      ar_of_arrays.first.each_with_index do |mz,i|
        if @mz_ranges.any? {|mz_range| mz_range === mz }
          new_mzs << mz
          new_int << ar_of_arrays.last[i]
        end
      end
      [new_mzs, new_int]
    else
      #puts "REJECTING: #{@rt_range} ===! #{time}"
      nil
    end
  end
end

filters = if opt[:filter].size == 0
            nil
          else
            opt[:filter].map {|string| Filter.from_string(string) } 
          end

# returns a scan number from an id string
def scan_num(id)
  md = /scan=(\d+)/.match(id)
  md[1].to_i if md
end

ARGV.each do |file|
  scan_nums = [] ; times = [] ; spectra = []
  MzML::Doc.open(file) do |mzml|
    mzml.index[:spectrum].keys.each do |id|
      spec = mzml.spectrum(id)
      if mz = spec.mz
        if filters
          ar_of_ars = filters.map {|filter| filter.filter(spec.retention_time.to_f,  [mz, spec.intensity])}.compact.first
          if ar_of_ars
            scan_nums << scan_num(spec.id)
            times << spec.retention_time.to_f
            spectra << ar_of_ars
          end
        else
          scan_nums << scan_num(spec.id)
          times << spec.retention_time.to_f
          spectra << [mz, spec.intensity]
        end
      end
    end
  end

  plms1 = Ms::Msrun::Plms1.new(scan_nums, times, spectra)
  plms1.write(file.chomp(File.extname(file)) + "#{filters && '.filtered'}" + '.plms1')
end
