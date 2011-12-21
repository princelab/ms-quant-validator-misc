#!/usr/bin/env ruby

require 'mzml'
require 'ms/msrun/plms1'
require 'yaml'
require 'trollop'

Point = Struct.new(:rt, :mz, :it, :rt_i, :mz_i)

class PIC
  attr_accessor :id
  attr_accessor :points
  attr_accessor :missing_points
  attr_accessor :out_of_range_points

  def initialize(id)
    @id = id
    (@points, @missing_points, @out_of_range_points) = [], [], []
  end
end

class FeatureMapper
  attr_accessor :rt_hash
  attr_accessor :mz_hashes

  attr_accessor :rt_range
  attr_accessor :mz_ranges

  def initialize(rt_hash, mz_hashes)
    (@rt_hash, @mz_hashes) = rt_hash, mz_hashes
    retention_times = rt_hash.keys.sort
    @rt_range = retention_times.min..retention_times.max
    @mz_ranges = @mz_hashes.map do |mz_hash| 
      sorted = mz_hash.keys.sort
      if (min=sorted.min) && (max=sorted.max)
        min..max
      end
    end
  end

  # returns a PIC file
  def create_pic(feature_mzml_file)
  end
end

  # takes a retention time, m/z and intensity triplet and converts to integer index
  def push(triplet)
    (rt, mz, intensity) = triplet
    if @rt_range === rt
      if @rt_hash.key?(rt)
        rt_index = @rt_hash[rt]
        mz_hash = @mz_hashes[rt_index]
        if mz_hash.key?(mz)
          @rt_indices << rt_index  # only push on if found m/z too
          @mz_indices << mz_hash[mz]
          @rts << rt
          @mzs << mz
          @intensities << intensity
        else
          if @mz_ranges[rt_index] === mz
            @rt_indices_missing << rt_index
            @mz_indices_missing << mz
          else
            @rt_indices_oor << rt_index
            @mz_indices_oor << mz
          end
        end
      else
        if @rt_range === rt
          @rt_indices_missing << rt
          @mz_indices_missing << mz
        else
          @rt_indices_oor << rt
          @mz_indices_oor << mz
        end
      end
    end
  end

  # returns the weighted average of the mzs
  def weighted_average
    mzs.zip(intensities).inject(0.0) {|sum,pair| sum += pair[0]*pair[1] } / intensities.reduce(:+)
  end

  # http://en.wikipedia.org/wiki/Mean_square_weighted_deviation
  # (the standard deviation is the square root of the variance)
  # returns the variacnce and the weighted avg
  def unbiased_weighted_estimate_of_variance
    # xbarstar is the weighted mean
    # ∑ is summation
    # ∑weights / ( (∑weights)^2 - ∑(weights^2) ) *?? ∑weight((x-xbarstar)^2)
    wavg = weighted_average
    sum_of_weights = intensities.reduce(:+)
    first_term = sum_of_weights.to_f / ( sum_of_weights**2 - intensities.map {|v| v**2 }.reduce(:+) ) 
    last_term = intensities.zip(mzs).map {|inten, mz|  inten * ( mz - wavg )**2 }.reduce(:+)
    [first_term * last_term, wavg]
  end

  def remove_outliers!
    (variance, wavg) = unbiased_weighted_estimate_of_variance
    stddev = Math.sqrt(variance)
    up_bound = wavg + (3*stddev)
    lo_bound = wavg - (3*stddev)
    keep = @rt_indices.zip(@mz_indices, @rts, @mzs, @intensities).map.select {|ar| ((ar[3] <= up_bound) && (ar[3] >= lo_bound)) }
    unfold_arrays!(keep)
  end

  # given a list of arrays structured: rt_index, mz_index, rt, mz, intensity
  # sets the appropriate instance variables.
  def unfold_arrays!(ar_of_data_ars)
    (@rt_indices, @mz_indices, @rts, @mzs, @intensities) = [], [], [], [], []
    instance_ars = [@rt_indices, @mz_indices, @rts, @mzs, @intensities]
    ar_of_data_ars.each do |ar|
      instance_ars.zip(ar) do |iv_ar, val|
        iv_ar << val
      end
    end
  end

  def remove_duplicates_in_scan!
    wavg = weighted_average 
    #               [0              1            2     3     4           ]
    rti_to_arrays = @rt_indices.zip(@mz_indices, @rts, @mzs, @intensities).group_by {|row| row[0] }
    keep = rti_to_arrays.map do |rt_i, data_ars|
      if data_ars.size == 1
        data_ars.first
      else
        data_ars.sort_by {|data_ar| (data_ar[3] - wavg).abs }.first
      end
    end
    unfold_arrays!(keep)
  end

  def inspect
    "<PICIndex #{rt_indices.zip(mz_indices).to_a.inspect}"
  end
end

parser = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} file.plms1 <feature>.mzML ...
output: pic_list.yml
}
  opt :remove_duplicates_in_scan, "if > 1 m/z per scan, retain the closest scan to the weighted average of all peaks", :default => false
  opt :remove_outliers, "remove anything > 3 standard deviations from the weighted mean", :default => false
end

opt = parser.parse(ARGV)

if ARGV.size < 2
  parser.educate && exit
end

plms1_file = ARGV.shift

rt_hash = {}  # float => integer index
plms1 = Ms::Msrun::Plms1.new.read(plms1_file)


File.open(plms1_file + 'a','w') do |out|
  plms1.times.zip(plms1.scan_numbers, plms1.spectra) do |time, scan_num, spectrum|
    out.puts "#{time}: #{scan_num}"
    out.puts spectrum.first.join(", ")
    out.puts spectrum.last.join(", ")
  end
end

plms1.times.each_with_index do |time,i|
  rt_hash[time] = i
end

# [ {float => mz index}, ... ]
mz_hashes = plms1.spectra.map do |spectrum|
  mz_hash = {}
  spectrum.first.each_with_index do |mz,i|
    mz_hash[mz] = i
  end
  mz_hash
end

pics = ARGV.map do |feature_file|
  pi = PICIndex.new(File.basename(feature_file,'.*'), rt_hash, mz_hashes)
  MzML::Doc.open(feature_file) do |mzml|
    mzml.index[:spectrum].keys.each do |id|
      spec = mzml.spectrum(id)
      mzs = spec.mz
      if mzs
        mzs.zip(spec.intensity) do |mz, inten|
          pi.push [spec.retention_time.to_f, mz, inten]
        end
      end
    end
  end
  pi
end

pics.select! {|pic| pic.rt_indices.size > 0 }

if opt.remove_duplicates_in_scan
  pics.each(&:remove_duplicates_in_scan!)
end
if opt.remove_outliers
  pics.each(&:remove_outliers!)
end

base = plms1_file.chomp(File.extname(plms1_file))
['', '_oor', '_missing'].each do |ext|
  hash = {}
  pics.each do |pic|
    doublets = pic.send("rt_indices#{ext}".to_sym).zip(pic.send("mz_indices#{ext}".to_sym)).map.to_a
    hash[pic.id] = doublets
  end

  File.open(base + ".featurePIC#{ext}.yml",'w') do |out|
    out.print hash.to_yaml
  end
  if ext == ''
    File.open(base +".featurePIC#{ext}.txt",'w') do |out| 
      hash.each do |id, pairs|
        out.puts "ID=#{id}"
        pairs.each do |pair|
          out.puts pair.join(" ")
        end
      end
    end
  end
end

if $0 == __FILE__
  require 'spec/more'
  describe 'tests' do

  end

end
